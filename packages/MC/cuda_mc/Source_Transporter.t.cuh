//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Transporter.t.cuh
 * \author Thomas M. Evans
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_Transporter_t_cuh
#define cuda_mc_Source_Transporter_t_cuh

#include <iomanip>
#include <iostream>
#include <cmath>
#include <chrono>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/count.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Launch_Args.t.cuh"
#include "cuda_utils/Work_Pool.cuh"
#include "harness/Diagnostics.hh"
#include "utils/String_Functions.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "Source_Transporter.hh"

#include "Particle_Vector.cuh"
#include "Physics.cuh"
#include "Domain_Transporter.cuh"
#include "VR_Roulette.cuh"
#include "Tallier.cuh"
#include "RNG_Control.cuh"
#include "Source_Provider.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// Functor to transport source particles
//---------------------------------------------------------------------------//
template <class Geometry>
class Transport_Functor
{
  public:

    typedef Domain_Transporter<Geometry>    Transporter_t;
    typedef Particle_Vector<Geometry>       Particle_Vector_t;

    // Constructor
    Transport_Functor( const Transporter_t     *trans,
                             Particle_Vector_t  particles,
                             cuda::Work_Pool   *pool )
        : d_transporter(trans)
        , d_particles(particles)
        , d_pool(pool)
    {
    }

    // Operator apply for functor (transport 1 particle)
    __device__ void operator()( std::size_t tid )
    {
        // Initialize work pool
        d_pool->initialize();

        while (d_pool->work_per_thread() > 0)
        {
            for (int i = 0; i < d_pool->work_per_thread(); ++i)
            {
                // Get particle index
                int pid = d_pool->work_id(i);
                if (pid >= 0 && pid < d_particles.size())
                {
                    DEVICE_CHECK(d_particles.alive(pid));

                    // transport the particle through this (replicated) domain
                    d_transporter->transport(pid,d_particles);

                    // Inactivate dead particles
                    if (!d_particles.alive(pid))
                        d_pool->set_inactive(i);
                }
            }

            // Consolidate
            d_pool->consolidate();
        }
    }

  private:

    // On-device pointers
    const Transporter_t *d_transporter;
    Particle_Vector_t    d_particles;
    cuda::Work_Pool     *d_pool;
};

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor/
 */
template <class Geometry>
Source_Transporter<Geometry>::Source_Transporter(RCP_Std_DB   db,
                                                 SDP_Geometry geometry,
                                                 SDP_Physics  physics)
    : d_geometry(geometry)
    , d_physics(physics)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
    , d_transport_time(0.0)
    , d_source_time(0.0)
{
    REQUIRE(!db.is_null());
    REQUIRE(d_geometry.get_host_ptr());
    REQUIRE(d_geometry.get_device_ptr());
    REQUIRE(d_physics.get_host_ptr());
    REQUIRE(d_physics.get_device_ptr());

    // Build variance reduction
    std::string var = profugus::to_lower(
        db->get<std::string>("variance reduction",std::string("roulette")) );
    if( var == "roulette" )
    {
        d_vr = cuda::shared_device_ptr<VR_Roulette_t>(db);
    }

    int seed = db->get("seed",1234);

    d_block_size = db->get("block_size",256);

    d_np_per_thread = db->get("np_per_thread",1);

    std::string verb = profugus::to_lower(db->get<std::string>("verbosity",
        std::string("none")));
    if (verb == "none")
        d_verbosity = NONE;
    else if (verb == "low")
        d_verbosity = LOW;
    else if (verb == "medium")
        d_verbosity = MEDIUM;
    else if (verb == "high")
        d_verbosity = HIGH;
    else
        INSIST(false,"Invalid verbosity.");

    d_rng_control = std::make_shared<RNG_Control>(seed);

    // Build domain transporter
    d_transporter = std::make_shared<Transporter_DMM_t>(
        db, d_geometry, d_physics, d_vr );

    // Build particle vector
    d_particle_vec = std::make_shared<Particle_Vector_DMM_t>();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set new tallier
 */
template <class Geometry>
void Source_Transporter<Geometry>::set(SP_Tallier_DMM tallier)
{
    REQUIRE(tallier);
    d_tallier = tallier;
    
    // Don't build on-device tallier here, wait until just before kernel launch
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 */
template <class Geometry>
void Source_Transporter<Geometry>::solve(SP_Source source) const
{
    // No global barriers in here.
    // Introducing sync points here will be detrimental to efficiency
    //  of particle batching
    REQUIRE(source);

    SCOPED_TIMER("MC::Source_Transporter.solve");

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> diff;

    // Build Tallier and Domain_Transporter
    auto sdp_tallier = cuda::shared_device_ptr<Tallier<Geometry>>(
        d_tallier->device_instance());
    d_transporter->set(sdp_tallier);
    auto sdp_transporter =
        cuda::shared_device_ptr<Domain_Transporter<Geometry>>(
            d_transporter->device_instance());

    // Get number of particles in source
    size_type num_particles      = source->num_batch();
    size_type num_particles_left = num_particles;

    if (d_verbosity >= LOW)
    {
        std::cout << "Starting solve with " << num_particles
            << " particles on node " << profugus::node() << std::endl;
    }

    thrust::device_vector<int> event(num_particles);
    thrust::device_vector<int> indirection(num_particles);
    thrust::counting_iterator<int> cnt(0);
    thrust::copy(cnt,cnt+num_particles,indirection.begin());

    // Get source particles
    start = std::chrono::high_resolution_clock::now();
    Source_Provider<Geometry> provider;
    provider.get_particles(source, d_rng_control, d_particle_vec, indirection);
    cudaDeviceSynchronize();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    d_source_time += diff.count();

    // Build launch args
    cuda::Launch_Args<cuda::arch::Device> launch_args;
    launch_args.set_block_size(d_block_size);

    // Compute number of threads from surviving particle count
    int num_threads = num_particles_left / d_np_per_thread;

    // Make launch a multiple of cuda warp size
    if (num_threads % 32)
        num_threads = (num_particles_left / 32 + 1) * 32;
    launch_args.set_num_elements(num_threads);

    // Build work pool
    int warps_per_block = launch_args.block_size() / 32;
    int num_warps = launch_args.grid_size() * warps_per_block;
    cuda::Work_Pool_DMM pool_dmm(indirection,num_warps);
    auto pool = cuda::shared_device_ptr<cuda::Work_Pool>(
        pool_dmm.device_instance());

    // Build and execute kernel
    cudaDeviceSynchronize();
    start = std::chrono::high_resolution_clock::now();
    Transport_Functor<Geometry> f(
            sdp_transporter.get_device_ptr(), 
            d_particle_vec->device_instance(),
            pool.get_device_ptr());

    cuda::parallel_launch( f, launch_args );
    REQUIRE(cudaSuccess == cudaGetLastError());
    cudaDeviceSynchronize();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    d_transport_time += diff.count();

    // increment the particle counter
    DIAGNOSTICS_ONE(integers["particles_transported"] += num_particles);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set a fission site container and keff for sampling fission sites.
 *
 * Setting a fission site container tells the fixed-source solver to sample
 * fission sites that can be used in an outer k-code calculation.  Thus, this
 * function should be called if the fixed-source solver is used as the inner
 * part of a k-code eigenvalue calculation.  It should be called once per
 * k-code iteration to update the eigenvalue.
 *
 * Fission sites are added to the container, it is \b not emptied.
 */
template <class Geometry>
void Source_Transporter<Geometry>::sample_fission_sites(
    SP_Fission_Site_Vec fis_sites, double keff)
{
    // set the transporter with the fission site container and the latest keff
    // iterate
    d_transporter->set(fis_sites, keff);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get number of fission sites sampled during last transport.
 *
 * The size of the fission site vector is not changed during transport
 * (dynamic reallocation inside a kernel is slow and dangerous), so we 
 * need to resize the vector to the number of sites that were actually created.
 * Here we get the number of sampled sites from the Domain_Transporter.
 * Note that this requires updating the host-side copy of the
 * Domain_Transporter because the number of fission sites is accumulated as
 * a member variable during transport.
 *
 */
template <class Geometry>
int Source_Transporter<Geometry>::num_sampled_fission_sites()
{
    return d_transporter->num_sampled_fission_sites();
}


} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.t.cuh
//---------------------------------------------------------------------------//
