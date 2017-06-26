//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.t.cuh
 * \author Thomas M. Evans
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_t_cuh
#define cuda_mc_Uniform_Source_t_cuh

#include <Teuchos_Array.hpp>

#include "harness/Soft_Equivalence.hh"

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Memory.cuh"
#include "cuda_utils/Hardware.hh"
#include "cuda_utils/Utility_Functions.hh"

#include "comm/Timing.hh"
#include "comm/global.hh"

#include "Uniform_Source.hh"

#include <algorithm>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Sample the source.
template<class Geometry, class Shape>
__global__
void sample_source_kernel( const Geometry* geometry,
			   const Shape* shape,
			   const int num_particle,
			   const double wt,
			   const double* erg_cdf,
			   const int num_group,
			   const int num_run,
			   const int batch_size,
			   Particle_Vector<Geometry>* particles )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int* indices = particles->event_indices( events::DEAD );

    if ( idx < num_particle )
    {
	// Get the particle index.
	int pidx = indices[idx];
    DEVICE_CHECK( particles->event(pidx) == events::DEAD );

	// sample the angle isotropically
	cuda_utils::Space_Vector omega;
	cuda_utils::utility::sample_angle( 
	    omega, particles->ran(pidx), particles->ran(pidx) );

	// sample the geometry shape to get a starting position
	cuda_utils::Space_Vector r = shape->sample(
	    particles->ran(pidx), particles->ran(pidx), particles->ran(pidx) );

	// intialize the geometry state
	geometry->initialize( r, omega, particles->geo_state(pidx) );

	// get the material id
	unsigned int matid = geometry->matid( particles->geo_state(pidx) );

	// initialize the physics state by manually sampling the group
	int group = cuda_utils::utility::sample_discrete_CDF(
	    num_group, erg_cdf, particles->ran(pidx) );
	DEVICE_CHECK( group < num_group );
	particles->set_group( pidx, group );

	// set the material id in the particle
	particles->set_matid( pidx, matid );

	// set particle weight
	particles->set_wt( pidx, wt );

	// set the event.
	particles->set_event( pidx, events::TAKE_STEP );

	// set the batch.
	int batch = (num_run + idx) / batch_size;
	particles->set_batch( pidx, batch );

    // sample the initial distance to collision
    particles->set_dist_mfp( pidx, -std::log(particles->ran(pidx)) );

	// make particle alive
	particles->live( pidx );
    }
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry, class Shape>
Uniform_Source<Geometry,Shape>::Uniform_Source(
    const RCP_Std_DB& db,
    const cuda_utils::Shared_Device_Ptr<Geometry>& geometry,
    const int num_groups,
    const int num_batch )
    : d_geometry( geometry )
    , d_num_groups( num_groups )
    , d_num_batch( num_batch )
    , d_np_requested(0)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(1.0)
    , d_np_left(0)
    , d_np_run(0)
{
    REQUIRE(!db.is_null());

    // store the total number of requested particles
    d_np_requested = 1000;
    if (db->isType<int>("Np"))
        d_np_requested = db->get<int>("Np");
    else if (db->isType<size_type>("Np"))
        d_np_requested = db->get<size_type>("Np");
    else if (db->isParameter("Np"))
        VALIDATE(false,"Unrecognized type for parameter Np.");

    // initialize the total
    d_np_total = d_np_requested;

    // get the spectral shape
    const auto &shape = db->get(
        "spectral_shape", 
	Teuchos::Array<double>(d_num_groups, 1.0));
    CHECK(shape.size() == d_num_groups);

    // calculate the normalization
    double norm = std::accumulate(shape.begin(), shape.end(), 0.0);
    CHECK(norm > 0.0);

    // assign to the shape cdf
    Teuchos::Array<double> host_cdf( d_num_groups, 0.0 );
    norm = 1.0 / norm;
    host_cdf[0] = shape[0] * norm;
    for ( int n = 1; n < host_cdf.size(); ++n )
    {
        host_cdf[n] = host_cdf[n-1] + shape[n] * norm;
    }
    CHECK(profugus::soft_equiv(1.0, host_cdf.back(), 1.0e-6));

    // Allocate and copy the cdf to the device.
    cuda_utils::memory::Malloc( d_erg_cdf, d_num_groups );
    cuda_utils::memory::Copy_To_Device( 
	d_erg_cdf, host_cdf.getRawPtr(), d_num_groups );

    // initialize timers in this class, which may be necessary because domains
    // with no source will not make this timer otherwise
#if UTILS_TIMING > 0
    profugus::Timing_Diagnostics::update_timer(
        "CUDA_MC::Uniform_Source.get_particles", 0.0);
#endif
}

//---------------------------------------------------------------------------//
/*
 * \brief Destructor.
 */
template <class Geometry, class Shape>
Uniform_Source<Geometry,Shape>::~Uniform_Source()
{
    cuda_utils::memory::Free( d_erg_cdf );
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source.
 * \param geometric_shape
 */
template <class Geometry, class Shape>
void Uniform_Source<Geometry,Shape>::build_source(
    const cuda_utils::Shared_Device_Ptr<Shape>& shape )
{
    REQUIRE(shape.get_host_ptr());

    SCOPED_TIMER("cuda_profugus::Uniform_Source.build_source");

    // store the spatial shape
    d_shape = shape;

    // build the source based on domain replication
    build_DR();

    // set counters
    d_np_left = d_np_domain;
    d_np_run  = 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a particle from the source.
 */
template <class Geometry, class Shape>
void Uniform_Source<Geometry,Shape>::get_particles(
    cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles )
{
    REQUIRE(d_wt > 0.0);
    REQUIRE(d_shape.get_device_ptr());

    // do nothing if no source
    if (!d_np_left)
    {
        return;
    }

    SCOPED_TIMER("CUDA_MC::Uniform_Source.get_particles");

    // Get the particles that are dead.
    size_type num_particle = particles.get_host_ptr()->get_event_size( events::DEAD );

    // Calculate the total number of particles we will create.
    size_type num_to_create = std::min( d_np_left, num_particle );

    if (num_to_create > 0)
    {
        // Get CUDA launch parameters.
        REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
        unsigned int threads_per_block = 
        cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
        unsigned int num_blocks = num_to_create / threads_per_block;
        if ( num_to_create % threads_per_block > 0 ) ++num_blocks;

        // Create the particles.
        sample_source_kernel<<<num_blocks,threads_per_block>>>(
            d_geometry.get_device_ptr(),
            d_shape.get_device_ptr(),
            num_to_create,
            d_wt,
            d_erg_cdf,
            d_num_groups,
            d_np_run,
            d_batch_size,
            particles.get_device_ptr() );

        cudaDeviceSynchronize();
        CHECK(cudaSuccess == cudaGetLastError());

        // update counters
        d_np_left -= num_to_create;
        d_np_run += num_to_create;
    }
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build domain replicated source.
 */
template <class Geometry, class Shape>
void Uniform_Source<Geometry,Shape>::build_DR()
{
    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / profugus::nodes();

    // increase the domain count to get an equivalent number per batch
    d_np_domain += d_np_domain % d_num_batch;
    CHECK( 0 == d_np_domain % d_num_batch );

    // Get the actual batch size.
    d_batch_size = d_np_domain / d_num_batch;

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total  = d_np_domain * profugus::nodes();
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Uniform_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.t.cuh
//---------------------------------------------------------------------------//
