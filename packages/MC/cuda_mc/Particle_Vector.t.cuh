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

#include "comm/Timing.hh"

#include <vector>

#include <thrust/device_ptr.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/remove.h>

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
    : d_original_size( num_particle )
    , d_size( num_particle )
    , d_matid_vec( d_size )
    , d_matid( d_matid_vec.data().get() )
    , d_group_vec( d_size )
    , d_group( d_group_vec.data().get() )      
    , d_wt_vec( d_size )
    , d_wt( d_wt_vec.data().get() )
    , d_rng_vec( d_size )
    , d_rng( d_rng_vec.data().get() )
    , d_alive_vec( d_size )
    , d_alive( d_alive_vec.data().get() )
    , d_geo_state_vec( d_size )
    , d_geo_state( d_geo_state_vec.data().get() )
    , d_event_vec( d_size )
    , d_event( d_event_vec.data().get() )
    , d_batch_vec( d_size )
    , d_batch( d_batch_vec.data().get() )
    , d_step_vec( d_size )
    , d_step( d_step_vec.data().get() )
    , d_dist_mfp_vec( d_size )
    , d_dist_mfp( d_dist_mfp_vec.data().get() )
    , d_event_indices_vec( events::END_EVENT*d_size )
    , d_event_indices( d_event_indices_vec.data().get() )
    , d_event_sizes( events::END_EVENT, 0 )
    , d_event_lid( d_size )
    , d_event_stencil( d_size )
    , d_event_steering( d_size )
{
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
    thrust::device_vector<int> device_seeds( host_seeds );

    // Initialize the generators.
    init_rng_kernel<<<num_blocks,threads_per_block>>>( 
	d_size, device_seeds.data().get(), d_rng );

    // Initialize all particles to DEAD.
    init_event_kernel<<<num_blocks,threads_per_block>>>( 
        d_size, d_event, d_event_indices );

    // All particles right now are dead.
    d_event_sizes[ events::DEAD ] = d_size;

    // Setup the local id vector for sorting.
    thrust::device_vector<int> lid( d_size, 1 );
    thrust::exclusive_scan( lid.begin(), lid.begin() + d_size, d_event_lid.begin() );

    // Do the first sort to initialize particle event state.
    sort_by_event( d_size );
}

//---------------------------------------------------------------------------//
// Sort the local indices by event key.
template <class Geometry>
void Particle_Vector<Geometry>::sort_by_event( const int sort_size )
{
    REQUIRE( sort_size <= d_size );

    // Create a stencil functor.
    Stencil_Functor stencil_functor;

    // Gather the indices for each event.
    for ( int e = 0; e < events::END_EVENT; ++e )
    {
        // Create the stencil vector.
        stencil_functor.d_event = static_cast<events::Event>(e);
        thrust::transform( thrust::device,
                           d_event_vec.begin(),
                           d_event_vec.begin()+d_size, 
                           d_event_stencil.begin(),
                           stencil_functor );

        // Get the number of particles with this event.
        d_event_sizes[e] = thrust::reduce( thrust::device,
                                           d_event_stencil.begin(),
                                           d_event_stencil.begin() + d_size );

        // If some particles had this event then extract their indices.
        if ( d_event_sizes[e] > 0 )
        {
            // Create the steering vector.
            thrust::exclusive_scan( thrust::device,
                                    d_event_stencil.begin(),
                                    d_event_stencil.begin() + d_size,
                                    d_event_steering.begin() );

            // Scatter the particles into the event indices. Indices will
            // always be strided by the original vector size to ensure
            // consistency between host and device.
            thrust::scatter_if(
                thrust::device,
                d_event_lid.begin(),
                d_event_lid.begin() + d_size,
                d_event_steering.begin(),
                d_event_stencil.begin(),
                d_event_indices_vec.begin() + e*d_original_size );
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
// Clear dead particles from the vector.
template<class Geometry>
void Particle_Vector<Geometry>::cultivate()
{
    // Create a stencil functor for dead particles.
    Stencil_Functor dead_stencil;
    dead_stencil.d_event = events::DEAD;

    // Cultivate the dead.
    thrust::remove_if( d_matid_vec.begin(), 
                       d_matid_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_group_vec.begin(), 
                       d_group_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_wt_vec.begin(), 
                       d_wt_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_rng_vec.begin(), 
                       d_rng_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_alive_vec.begin(), 
                       d_alive_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_geo_state_vec.begin(), 
                       d_geo_state_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_batch_vec.begin(), 
                       d_batch_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_step_vec.begin(), 
                       d_step_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    thrust::remove_if( d_dist_mfp_vec.begin(), 
                       d_dist_mfp_vec.begin() + d_size, 
                       d_event_vec.begin(), 
                       dead_stencil );

    // Events are last so we don't destroy the stencil until were are done
    // with it.
    auto event_it = thrust::remove_if( d_event_vec.begin(),
                                       d_event_vec.begin() + d_size,
                                       dead_stencil );
    
    // Update the current vector size.
    d_size = thrust::distance( d_event_vec.begin(), event_it );
}

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // end cuda_mc_Particle_Vector_t_cuh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.t.cuh
//---------------------------------------------------------------------------//
