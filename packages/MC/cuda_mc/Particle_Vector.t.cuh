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
#include "cuda_utils/Utility_Functions.hh"

#include "comm/Timing.hh"

#include <vector>

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
// Initialize particle local ids.
__global__ void init_lid_kernel( const int size, int* lids )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) lids[idx] = idx;
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
// Get event lower bounds and clear the number of events.
__global__ void event_bounds_kernel( const int size,
                                     int* num_event,
                                     int* event_bounds )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( idx < size )
    {
        // Get the lower bound.
        event_bounds[idx] = 0;
        for ( int i = 0; i < idx; ++i ) 
        {
            event_bounds[idx] += num_event[i];
        }

        // Clear the number of events.
        num_event[idx] = 0;
    }
}

//---------------------------------------------------------------------------//
// Reorder the local ids.
__global__ void reorder_lid_kernel( const int sort_size,
                                    const int vector_size,
                                    const int* event_bounds,
                                    const int* event_bins,
                                    const events::Event* events,
                                    int* lids )
{
  DEVICE_REQUIRE( sort_size <= vector_size );

  // Get the thread index.
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if ( idx < sort_size )
  {
      // Get the event corresponding to this vector index.
      const int* event_lb = cuda_utils::utility::lower_bound(
          event_bounds,
          event_bounds + events::END_EVENT,
          idx,
          cuda_utils::utility::upper_bound_comp<int>() );

      int event = event_lb - event_bounds - 1;
      DEVICE_CHECK( event < events::END_EVENT );
      DEVICE_CHECK( idx >= event_bounds[event] );

      int local_bin_index = idx - event_bounds[event];
      DEVICE_CHECK( local_bin_index < vector_size );

      // Copy the particle index from the bin into the local index array.      
      lids[idx] = event_bins[ event*vector_size + local_bin_index ];
      DEVICE_CHECK( event == events::DEAD || event == events[lids[idx]] );
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
    cuda_utils::memory::Malloc( d_matid, d_size );
    cuda_utils::memory::Malloc( d_group, d_size );
    cuda_utils::memory::Malloc( d_wt, d_size );
    cuda_utils::memory::Malloc( d_rng, d_size );
    cuda_utils::memory::Malloc( d_alive, d_size );
    cuda_utils::memory::Malloc( d_geo_state, d_size );
    cuda_utils::memory::Malloc( d_event, d_size );
    cuda_utils::memory::Malloc( d_lid, d_size );
    cuda_utils::memory::Malloc( d_batch, d_size );
    cuda_utils::memory::Malloc( d_step, d_size );
    cuda_utils::memory::Malloc( d_dist_mfp, d_size );
    cuda_utils::memory::Malloc( d_event_bounds, events::END_EVENT );
    cuda_utils::memory::Malloc( d_num_event, events::END_EVENT );
    cuda_utils::memory::Malloc( d_event_bins, events::END_EVENT*d_size );    

    // Get CUDA launch parameters.
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    // Create seeds for the random number generators.
    std::vector<int> host_seeds( d_size );
    for ( auto& s : host_seeds ) s = rng.uniform<int>();

    // Copy the seeds to the device.
    int* device_seeds = NULL;
    cuda_utils::memory::Malloc( device_seeds, d_size );
    cuda_utils::memory::Copy_To_Device( device_seeds, host_seeds.data(), d_size );

    // Initialize the generators.
    init_rng_kernel<<<num_blocks,threads_per_block>>>( 
	d_size, device_seeds, d_rng );

    // Deallocate the device seeds.
    cuda_utils::memory::Free( device_seeds );

    // Initialize the number of events.
    clear_num_event_kernel<<<1,events::END_EVENT>>>(
        events::END_EVENT, d_num_event );

    // Initialize all particles to DEAD.
    init_event_kernel<<<num_blocks,threads_per_block>>>( 
        d_size, d_event, d_num_event, d_event_bins );
    d_event_sizes[ events::DEAD ] = d_size;
    
    // Do the first sort to initialize particle event state.
    sort_by_event( d_size );
}
    
//---------------------------------------------------------------------------//
// Destructor.
template <class Geometry>
Particle_Vector<Geometry>::~Particle_Vector()
{
    cuda_utils::memory::Free( d_matid );
    cuda_utils::memory::Free( d_group );
    cuda_utils::memory::Free( d_wt );
    cuda_utils::memory::Free( d_rng );
    cuda_utils::memory::Free( d_alive );
    cuda_utils::memory::Free( d_geo_state );
    cuda_utils::memory::Free( d_event );
    cuda_utils::memory::Free( d_lid );
    cuda_utils::memory::Free( d_batch );
    cuda_utils::memory::Free( d_step );
    cuda_utils::memory::Free( d_dist_mfp );
    cuda_utils::memory::Free( d_event_bounds );
    cuda_utils::memory::Free( d_num_event );
    cuda_utils::memory::Free( d_event_bins );
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
    SCOPED_TIMER("CUDA_MC::Particle_Vector.sort_by_event");

    REQUIRE( sort_size <= d_size );

    // Get CUDA launch parameters.
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
    unsigned int num_blocks = sort_size / threads_per_block;
    if ( sort_size % threads_per_block > 0 ) ++num_blocks;

    // Get the number of particles with each event.
    cuda_utils::memory::Copy_To_Host( d_event_sizes.getRawPtr(),
                                      d_num_event,
                                      events::END_EVENT );

    // Get the event bounds.
    unsigned int num_bounds_blocks = events::END_EVENT / threads_per_block;
    if ( events::END_EVENT % threads_per_block > 0 ) ++num_bounds_blocks;
    event_bounds_kernel<<<num_bounds_blocks,threads_per_block>>>(
        static_cast<int>(events::END_EVENT), d_num_event, d_event_bounds );

    if (sort_size > 0)
    {
        // Reorder the local ids.
        reorder_lid_kernel<<<num_blocks,threads_per_block>>>( sort_size,
                                                              d_size,
                                                              d_event_bounds,
                                                              d_event_bins,
                                                              d_event,
                                                              d_lid );
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

} // end namespace profugus

#endif // end cuda_mc_Particle_Vector_t_cuh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.t.cuh
//---------------------------------------------------------------------------//
