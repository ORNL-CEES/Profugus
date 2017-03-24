//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Keff_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Wed May 14 13:29:40 2014
 * \brief  Keff_Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Keff_Tally_t_cuh
#define cuda_mc_Keff_Tally_t_cuh

#include "Keff_Tally.hh"

#include "comm/global.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Memory.cuh"
#include "cuda_utils/Hardware.hh"
#include "Definitions.hh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Reset keff device
__global__ void reset_keff_kernel( const int num_particles,
				   double* keff )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( idx < num_particles )
    {
	keff[idx] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Tally accumulate kernel.
template<class Geometry>
__global__ void accumulate_kernel( const Physics<Geometry>* phyiscs,
				   const Particle_Vector<Geometry>* particles,
				   const std::size_t num_collision,
				   const std::size_t num_boundary,
				   double* keff )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int* collision_indices = particles->event_indices( events::COLLISION );
    int* boundary_indices = particles->event_indices( events::BOUNDARY );

    if ( idx < num_collision + num_boundary )
    {
	// Get the particle index.
	int pidx = ( idx < num_collision )
		   ? collision_indices[idx]
		   : boundary_indices[idx - num_collision];
    
	// Tally keff
	keff[idx] += particles->wt(pidx) * particles->step(pidx) *
                     phyiscs->total( physics::NU_FISSION,
                                     particles->matid(pidx),
                                     particles->group(pidx) );
    }
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Kcode solver should construct this with initial keff estimate.
 */
template <class Geometry>
Keff_Tally<Geometry>::Keff_Tally(
    const double keff_init,
    const cuda_utils::Shared_Device_Ptr<Physics_t>& physics,
    const int vector_size )
    : d_physics( physics )
    , d_keff_cycle(keff_init)
    , d_keff_device(nullptr)
    , d_vector_size( vector_size )
{
    REQUIRE(physics);

    // Allocate the keff work vectors.
    d_keff_host.resize( d_vector_size );
    cuda_utils::memory::Malloc( d_keff_device, d_vector_size );
    
    // reset tally
    reset();

    // Reset device cycle keff
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
        cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
    unsigned int num_blocks = d_vector_size / threads_per_block;
    if ( d_vector_size % threads_per_block > 0 ) ++num_blocks;
    reset_keff_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>(
        d_vector_size,
        d_keff_device );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Geometry>
Keff_Tally<Geometry>::~Keff_Tally()
{
    cuda_utils::memory::Free( d_keff_device );
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate average keff over active cycles.
 *
 * This should really only be called after begin_active_cycles; but either way
 * it will return the average keff since either initialization of this tally or
 * the last call to begin_active_cycles.
 *
 * This is only meaningful after one complete cycle. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::mean() const
{
    if (d_cycle < 1)
        return -1.;

    double keff = d_keff_sum / d_cycle;

    ENSURE(keff >= 0.);
    return keff;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate variance of average keff estimate over active cycles.
 *
 * This is only meaningful after two complete cycles. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::variance() const
{
    if (d_cycle < 2)
        return d_keff_sum * d_keff_sum;

    double var = (d_keff_sum_sq - d_keff_sum * d_keff_sum / d_cycle)
                    / (d_cycle * (d_cycle - 1));

    ENSURE(var >= 0.);
    return var;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Track particle and do tallying.
 *
 * This only uses the particle's weight. The accumulated tally is
 * \f[
   k_{\mbox{eff}} = wl\nu\sigma_{\mbox{f}}
 * \f]
 * where \f$l\f$ is the step-length.
 */
template <class Geometry>
void Keff_Tally<Geometry>::accumulate(
    const cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    // Get the particles that just had a collision.
    int num_collision = 
        particles.get_host_ptr()->get_event_size( events::COLLISION );

    // Get the particles that just had a boundary.
    int num_boundary = 
        particles.get_host_ptr()->get_event_size( events::BOUNDARY );

    // Calculate the launch parameters.
    std::size_t num_particle = num_collision + num_boundary;
    if (num_particle > 0)
    {
        REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
        unsigned int threads_per_block = 
        cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
        unsigned int num_blocks = num_particle / threads_per_block;
        if ( num_particle % threads_per_block > 0 ) ++num_blocks;

        // Process the tallies.
        accumulate_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>(
        d_physics.get_device_ptr(),
        particles.get_device_ptr(),
        num_collision,
        num_boundary,
        d_keff_device );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin active cycles in a kcode calculation.
 *
 * This resets the accumulated keff statistics.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_active_cycles()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin a new keff cycle
 *
 * This clears the current accumulated path lengths.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_cycle()
{
    // Reset cycle keff
    d_keff_cycle = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Mark the end of a keff cycle
 *
 * This performs a global sum across processors, because the provided
 * num_particles is the total number over all domains (blocks plus sets).
 *
 * We accumulate the sum and sum-of-squares of the keff so that we can calculate
 * averages and variances.
 */
template <class Geometry>
void Keff_Tally<Geometry>::end_cycle(double num_particles)
{
    REQUIRE(num_particles > 0.);

    // Pull tallies off the device.
    cuda_utils::memory::Copy_To_Host_Async( 
	d_keff_host.data(), d_keff_device, d_vector_size, d_stream );

    // Synchronize on this thread.
    d_stream.synchronize();
    REQUIRE(cudaSuccess == cudaGetLastError());

    // Add them to cycle tally.
    for ( int n = 0; n < d_vector_size; ++n )
    {
	d_keff_cycle += d_keff_host[n];
    }

    // Keff estimate is total nu-sigma-f reaction rate / num particles
    d_keff_cycle /= num_particles;

    // Do a global sum (since num_particles is global)
    profugus::global_sum(d_keff_cycle);

    // Accumulate first and second moments of cycles, as well as counter
    ++d_cycle;
    d_keff_sum    += d_keff_cycle;
    d_keff_sum_sq += d_keff_cycle * d_keff_cycle;

    // Store keff estimate
    d_all_keff.push_back(d_keff_cycle);

    // Reset device cycle keff
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
        cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
    unsigned int num_blocks = d_vector_size / threads_per_block;
    if ( d_vector_size % threads_per_block > 0 ) ++num_blocks;

    reset_keff_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>(
        d_vector_size,
        d_keff_device );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clear accumulated keff values
 *
 * This will still preserve the latest keff value so that existing fission sites
 * will be sampled correctly.
 */
template <class Geometry>
void Keff_Tally<Geometry>::reset()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;

    d_all_keff.clear();
}

} // end namespace cuda_profugus

#endif // cuda_mc_Keff_Tally_t_cuh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.t.hh
//---------------------------------------------------------------------------//
