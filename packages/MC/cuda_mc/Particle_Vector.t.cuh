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

#include "comm/Timing.hh"

#include <vector>

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/distance.h>
#include <thrust/binary_search.h>

#include <cuda_runtime.h>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Initialize particle random number generators
__global__ void init_rng_kernel( const std::size_t size,
				 const int* seeds,
				 curandState* rng )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) curand_init( seeds[idx], 0, 0, &rng[idx] );
}

//---------------------------------------------------------------------------//
// Initialize particle local ids.
__global__ void init_lid_kernel( const std::size_t size, std::size_t* lids )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) lids[idx] = idx;
}

//---------------------------------------------------------------------------//
// Initialize particles to DEAD.
__global__ void init_event_kernel( const std::size_t size,
				   events::Event* events )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) events[idx] = events::DEAD;
}

//---------------------------------------------------------------------------//
// Reset event bounds.
__global__ void reset_event_bounds_kernel( int* event_bounds )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Initialize bounds.
    if ( idx < events::END_EVENT )
    {
        event_bounds[ 2*idx ] = 0;
        event_bounds[ 2*idx + 1 ] = 0;
    }
}

//---------------------------------------------------------------------------//
// Get event bounds.
__global__ void event_bounds_kernel( const std::size_t size,
                                     const events::Event* events,
                                     int* event_bounds )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Get the lower bound.
    if ( 0 == idx )
    {
        event_bounds[ 2*events[idx] ] = idx;
    }
    else if ( idx < size )
    {
        if ( events[idx-1] < events[idx] )
        {
            event_bounds[ 2*events[idx] ] = idx;
        }
    }
    
    // Get the upper bound.
    if ( idx == size - 1 )
    {
        event_bounds[ 2*events[idx] + 1 ] = idx + 1;
    }
    else if ( idx < size-1 )
    {
        if ( events[idx+1] > events[idx] )
        {
            event_bounds[ 2*events[idx] + 1] = idx + 1;
        }
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
    cuda::memory::Malloc( d_lid, d_size );
    cuda::memory::Malloc( d_batch, d_size );
    cuda::memory::Malloc( d_step, d_size );
    cuda::memory::Malloc( d_dist_mfp, d_size );
    cuda::memory::Malloc( d_event_bounds, 2*events::END_EVENT );

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

    // Create the local ids.
    init_lid_kernel<<<num_blocks,threads_per_block>>>( d_size, d_lid );

    // Initialize all particles to DEAD.
    init_event_kernel<<<num_blocks,threads_per_block>>>( d_size, d_event );

    // All particles right now are dead.
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
    cuda::memory::Free( d_lid );
    cuda::memory::Free( d_batch );
    cuda::memory::Free( d_step );
    cuda::memory::Free( d_dist_mfp );
    cuda::memory::Free( d_event_bounds );
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
    SCOPED_TIMER("CUDA_MC::Particle_Vector.sort_by_event");    
    REQUIRE( sort_size <= d_size );

    // Sort the events.
    thrust::device_ptr<Event_t> event_begin( d_event );
    thrust::device_ptr<Event_t> event_end( d_event + sort_size );
    thrust::device_ptr<std::size_t> lid_begin( d_lid );
    thrust::sort_by_key( thrust::device, event_begin, event_end, lid_begin );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    // Reset event bounds
    reset_event_bounds_kernel<<<1,events::END_EVENT>>>( d_event_bounds );

    // Bin them.
    event_bounds_kernel<<<num_blocks,threads_per_block>>>(
        d_size, d_event, d_event_bounds );

    // Calculate event sizes.
    int work_size = 2 * events::END_EVENT;
    Teuchos::Array<int> work( work_size );
    cuda::memory::Copy_To_Host( work.getRawPtr(),
                                d_event_bounds,
                                work_size );
    for ( int i = 0; i < events::END_EVENT; ++i )
    {
        d_event_sizes[ i ] = work[ 2*i + 1 ] - work[ 2*i ];
    }
}

//---------------------------------------------------------------------------//
// Get the number of particles with a given event on the host.
template<class Geometry>
int Particle_Vector<Geometry>::get_event_size( const events::Event event ) const
{
    return d_event_sizes[ event ];
}

//---------------------------------------------------------------------------//
// Reset the vector.
template<class Geometry>
void Particle_Vector<Geometry>::reset()
{
    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    // Initialize all particles to DEAD.
    init_event_kernel<<<num_blocks,threads_per_block>>>( d_size, d_event );
    d_event_sizes[ events::DEAD ] = d_size;
}

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // end cuda_mc_Particle_Vector_t_cuh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.t.cuh
//---------------------------------------------------------------------------//
