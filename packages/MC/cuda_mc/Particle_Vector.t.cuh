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
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

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
// Initialize particles to DEAD.
__global__ void init_event_kernel( const std::size_t size,
				   events::Event* events,
                                   int* event_indices )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) 
    {
        events[idx] = events::DEAD;
        event_indices[events::DEAD*size + idx] = idx;
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
    , d_event_lid( d_size )
    , d_event_stencil( d_size )
    , d_event_steering( d_size )
{
    // Allocate data arrays.
    cuda::memory::Malloc( d_matid, d_size );
    cuda::memory::Malloc( d_group, d_size );
    cuda::memory::Malloc( d_wt, d_size );
    cuda::memory::Malloc( d_rng, d_size );
    cuda::memory::Malloc( d_alive, d_size );
    cuda::memory::Malloc( d_geo_state, d_size );
    cuda::memory::Malloc( d_event, d_size );
    cuda::memory::Malloc( d_batch, d_size );
    cuda::memory::Malloc( d_step, d_size );
    cuda::memory::Malloc( d_dist_mfp, d_size );
    cuda::memory::Malloc( d_event_indices, events::END_EVENT*d_size );

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

    // Initialize all particles to DEAD.
    init_event_kernel<<<num_blocks,threads_per_block>>>( 
        d_size, d_event, d_event_indices );

    // All particles right now are dead.
    d_event_sizes[ events::DEAD ] = d_size;

    // Setup the local id vector for sorting.
    thrust::device_vector<int> lid( d_size, 1 );
    thrust::exclusive_scan( lid.begin(), lid.end(), d_event_lid.begin() );

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
    cuda::memory::Free( d_batch );
    cuda::memory::Free( d_step );
    cuda::memory::Free( d_dist_mfp );
    cuda::memory::Free( d_event_indices );
}

//---------------------------------------------------------------------------//
// Sort the local indices by event key.
template <class Geometry>
void Particle_Vector<Geometry>::sort_by_event( const int sort_size )
{
    REQUIRE( sort_size <= d_size );

    // Get pointers to the event array.
    thrust::device_ptr<events::Event> event_begin( d_event );
    thrust::device_ptr<events::Event> event_end( d_event + d_size );
    
    // Create a stencil functor.
    Stencil_Functor stencil_functor;

    // Gather the indices for each event.
    for ( int e = 0; e < events::END_EVENT; ++e )
    {
        // Create the stencil vector.
        stencil_functor.d_event = static_cast<events::Event>(e);
        thrust::transform( thrust::device,
                           event_begin, 
                           event_end, 
                           d_event_stencil.begin(),
                           stencil_functor );

        // Get the number of particles with this event.
        d_event_sizes[e] = thrust::reduce( thrust::device,
                                           d_event_stencil.begin(),
                                           d_event_stencil.end() );

        // If some particles had this event then extract their indices.
        if ( d_event_sizes[e] > 0 )
        {
            // Create the steering vector.
            thrust::exclusive_scan( thrust::device,
                                    d_event_stencil.begin(),
                                    d_event_stencil.end(),
                                    d_event_steering.begin() );

            // Scatter the particles into the event indices.
            thrust::scatter_if(
                thrust::device,
                d_event_lid.begin(),
                d_event_lid.end(),
                d_event_steering.begin(),
                d_event_stencil.begin(),
                thrust::device_ptr<int>(d_event_indices + e*d_size) );
        }
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
