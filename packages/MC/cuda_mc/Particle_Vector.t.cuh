//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle_Vector.t.cuh
 * \author Stuart Slattery
 * \brief  Particle vector class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Particle_Vector.hh"

#include "cuda_utils/Hardware.hh"

#include <vector>

#include <thrust/device_ptr.h>
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
__global__ void init_rng( const std::size_t start_idx,
			  const int* seeds,
			  curandState* rng )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x + start_idx;
    curand_init( seeds[idx], 0, 0, &rng[idx] );
}

//---------------------------------------------------------------------------//
// Initialize particle local ids.
__global__ void init_lid( const std::size_t start_idx,
			  std::size_t* lids )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x + start_idx;
    lids[idx] = idx;
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
{
    // Allocate data arrays.
    cudaMalloc( (void**) &d_matid, d_size * sizeof(int) );
    cudaMalloc( (void**) &d_group, d_size * sizeof(int) );
    cudaMalloc( (void**) &d_wt, d_size * sizeof(double) );
    cudaMalloc( (void**) &d_rng, d_size * sizeof(curandState) );
    cudaMalloc( (void**) &d_alive, d_size * sizeof(bool) );
    cudaMalloc( (void**) &d_geo_state, d_size * sizeof(Geo_State_t) );
    cudaMalloc( (void**) &d_event, d_size * sizeof(Event_t) );
    cudaMalloc( (void**) &d_lid, d_size * sizeof(std::size_t) );
    cudaMalloc( (void**) &d_batch, d_size * sizeof(int) );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    unsigned int remainder = d_size % threads_per_block;

    // Create seeds for the random number generators.
    std::vector<int> host_seeds( d_size );
    for ( auto& s : host_seeds ) s = rng.uniform<int>();

    // Copy the seeds to the device.
    int* device_seeds;
    cudaMalloc( (void**) &device_seeds, d_size * sizeof(int) );
    cudaMemcpy( device_seeds, host_seeds.data(), d_size * sizeof(int),
		cudaMemcpyHostToDevice );

    // Initialize the generators.
    if ( num_blocks > 0 )
    {
	init_rng<<<num_blocks,threads_per_block>>>( 0, device_seeds, d_rng );
    }
    if ( remainder > 0 )
    {
	init_rng<<<1,remainder>>>( 
	    threads_per_block * num_blocks, device_seeds, d_rng );
    }

    // Deallocate the device seeds.
    cudaFree( device_seeds );

    // Create the local ids.
    if ( num_blocks > 0 )
    {
    	init_lid<<<num_blocks,threads_per_block>>>( 0, d_lid );
    }
    if ( remainder > 0 )
    {
    	init_lid<<<1,remainder>>>( threads_per_block * num_blocks, d_lid );
    }
}
    
//---------------------------------------------------------------------------//
// Destructor.
template <class Geometry>
Particle_Vector<Geometry>::~Particle_Vector()
{
    cudaFree( d_matid );
    cudaFree( d_group );
    cudaFree( d_wt );
    cudaFree( d_rng );
    cudaFree( d_alive );
    cudaFree( d_geo_state );
    cudaFree( d_event );
    cudaFree( d_lid );
    cudaFree( d_batch );
}

//---------------------------------------------------------------------------//
// Sort the local indices by event key.
template <class Geometry>
void Particle_Vector<Geometry>::sort_by_event()
{
    thrust::device_ptr<Event_t> event_begin( d_event );
    thrust::device_ptr<Event_t> event_end( d_event + d_size );
    thrust::device_ptr<std::size_t> lid_begin( d_lid );
    thrust::sort_by_key( event_begin, event_end, lid_begin );
}

//---------------------------------------------------------------------------//
// Given an event, get the index at which it starts and the number of
// particles with that event.
template <class Geometry>
void Particle_Vector<Geometry>::get_event_particles( 
    const Event_t event, 
    std::size_t& start_index, 
    std::size_t& num_particle ) const
{
    thrust::device_ptr<const Event_t> event_begin( d_event );
    thrust::device_ptr<const Event_t> event_end( d_event + d_size );
    auto event_range = thrust::equal_range( event_begin, event_end, event );
    start_index = thrust::distance( event_begin, event_range.first );
    num_particle = thrust::distance( event_range.first, event_range.second );
}

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.t.cuh
//---------------------------------------------------------------------------//
