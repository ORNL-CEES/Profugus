//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle_Vector.t.cuh
 * \author Stuart Slattery
 * \brief  Particle vector class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Particle_Vector_t_cuh
#define cuda_mc_Particle_Vector_t_cuh

#include "Particle_Vector.hh"

#include "cuda_utils/Hardware.hh"
#include "cuda_utils/Memory.cuh"

#include <vector>

#include <thrust/device_ptr.h>
#include <thrust/copy.h>

#include <cuda_runtime.h>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Initialize particle random number generators
__global__ void init_rng_kernel( const int size,
				 const int* seeds,
				 curandState* rng )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) curand_init( seeds[idx], 0, 0, &rng[idx] );
}

//---------------------------------------------------------------------------//
// Initialize particles to DEAD.
__global__ void init_event_kernel( const int size,
				   events::Event* events,
                                   int* num_event,
                                   int* event_bins )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) 
    {
        events[idx] = events::DEAD;
        atomicAdd( &num_event[events::DEAD], 1 );
        event_bins[ events::DEAD*size + idx ] = idx;
    }
}

//---------------------------------------------------------------------------//
// Reset the first event indicator for the particles.
__global__ void reset_first_event_kernel( const int size,
                                          bool* first_event )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size )
    {
        first_event[idx] = true;
    }
}

//---------------------------------------------------------------------------//
// Clear the number of events on the device.
__global__ void clear_num_event_kernel( const int size,
                                         int* num_event )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size )
    {
        num_event[idx] = 0;
    }
}

//---------------------------------------------------------------------------//
// HOST API
//---------------------------------------------------------------------------//
/*
 * \brief Constructor
 */
template <class Geometry>
Particle_Vector<Geometry>::Particle_Vector( const int num_particle, 
					    const profugus::RNG& rng )
    : d_size( num_particle )
    , d_event_sizes( events::END_EVENT, 0 )
{
    // Allocate data arrays.
    cuda::memory::Malloc( d_matid, d_size );
    cuda::memory::Malloc( d_group, d_size );
    cuda::memory::Malloc( d_wt, d_size );
    cuda::memory::Malloc( d_rng, d_size );
    cuda::memory::Malloc( d_alive, d_size );
    cuda::memory::Malloc( d_geo_state, d_size );
    cuda::memory::Malloc( d_event, d_size );
    cuda::memory::Malloc( d_first_event, d_size );
    cuda::memory::Malloc( d_batch, d_size );
    cuda::memory::Malloc( d_step, d_size );
    cuda::memory::Malloc( d_dist_mfp, d_size );
    cuda::memory::Malloc( d_num_event, events::END_EVENT );
    cuda::memory::Malloc( d_last_event_bins, events::END_EVENT*d_size );
    cuda::memory::Malloc( d_current_event_bins, events::END_EVENT*d_size );        

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    // Create seeds for the random number generators.
    std::vector<int> host_seeds( d_size );
    for ( auto& s : host_seeds ) s = rng.uniform<int>();

    // Copy the seeds to the device.
    int* device_seeds = NULL;
    cuda::memory::Malloc( device_seeds, d_size );
    cuda::memory::Copy_To_Device( device_seeds, host_seeds.data(), d_size );

    // Initialize the generators.
    init_rng_kernel<<<num_blocks,threads_per_block>>>( 
	d_size, device_seeds, d_rng );

    // Deallocate the device seeds.
    cuda::memory::Free( device_seeds );

    // Initialize the number of events.
    clear_num_event_kernel<<<1,events::END_EVENT>>>(
        events::END_EVENT, d_num_event );

    // Reset the first event indicator.
    reset_first_event_kernel<<<num_blocks,threads_per_block>>>( 
        d_size, d_first_event );

    // Initialize all particles to DEAD.
    init_event_kernel<<<num_blocks,threads_per_block>>>( 
        d_size, d_event, d_num_event, d_current_event_bins );
    d_event_sizes[ events::DEAD ] = d_size;
    
    // Do the first sort to initialize particle event state.
    sort_by_event( d_size );
}
    
//---------------------------------------------------------------------------//
// Destructor.
template <class Geometry>
Particle_Vector<Geometry>::~Particle_Vector()
{
    cuda::memory::Free( d_matid );
    cuda::memory::Free( d_group );
    cuda::memory::Free( d_wt );
    cuda::memory::Free( d_rng );
    cuda::memory::Free( d_alive );
    cuda::memory::Free( d_geo_state );
    cuda::memory::Free( d_event );
    cuda::memory::Free( d_first_event );
    cuda::memory::Free( d_batch );
    cuda::memory::Free( d_step );
    cuda::memory::Free( d_dist_mfp );
    cuda::memory::Free( d_num_event );
    cuda::memory::Free( d_last_event_bins );
    cuda::memory::Free( d_current_event_bins );
}

//---------------------------------------------------------------------------//
// Return if the vector is empty.
template <class Geometry>
bool Particle_Vector<Geometry>::empty() const
{
    int num_left = 0;
    for ( auto e : d_event_sizes ) num_left += e;
    return (0 == num_left);
}

//---------------------------------------------------------------------------//
// Get the number of particles that are not dead.
template <class Geometry>
int Particle_Vector<Geometry>::num_alive() const
{
    int num_alive = 0;
    for ( auto i = 0; i < events::DEAD; ++i ) num_alive += d_event_sizes[i];
    return num_alive;
}

//---------------------------------------------------------------------------//
// Sort the local indices by event key.
template <class Geometry>
void Particle_Vector<Geometry>::sort_by_event( const int sort_size )
{
    REQUIRE( sort_size <= d_size );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    // Reset the first event indicator so every particle starts with a first
    // event next cycle.
    reset_first_event_kernel<<<num_blocks,threads_per_block>>>( 
        d_size, d_first_event );

    // Get the number of particles with each event.
    cuda::memory::Copy_To_Host( d_event_sizes.getRawPtr(),
                                d_num_event,
                                events::END_EVENT );

    // Clear the count for the next round.
    clear_num_event_kernel<<<1,events::END_EVENT>>>(
        events::END_EVENT, d_num_event );

    // Copy the event bins.
    cuda::memory::Copy_Device_To_Device( d_last_event_bins, 
                                         d_current_event_bins,
                                         d_size*events::END_EVENT );
}

//---------------------------------------------------------------------------//
// Get the number of particles with a given event on the host.
template<class Geometry>
int Particle_Vector<Geometry>::get_event_size( const events::Event event ) const
{
    return d_event_sizes[ event ];
}

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // end cuda_mc_Particle_Vector_t_cuh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.t.cuh
//---------------------------------------------------------------------------//
