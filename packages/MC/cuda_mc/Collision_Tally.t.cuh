//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Collision_Tally.t.cu
 * \author Stuart Slattery
 * \brief  Collision class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Collision_Tally_t_cuh
#define cuda_mc_Collision_Tally_t_cuh

#include "cuda_utils/CudaDBC.hh"

#include "mc/Definitions.hh"

#include "Collision_Tally.hh"

#include <cuda_runtime.h>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Initialize the tally to zero
__global__ void init_tally_kernel( const std::size_t size, std::size_t* tally )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) tally[idx] = 0.0;
}

//---------------------------------------------------------------------------//
// Tally particles that have had a collision.
__global__ void tally_kernel( const Geometry* geometry,
			      const Particle_Vector* particles,
			      const int start_idx,
			      const int num_particle,
			      const int num_batch,
			      const int num_cell,
			      double* tally )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < num_particle )
    {
	std::size_t pidx = idx + start_idx;
	std::size_t tally_idx = particles->batch( pidx ) * num_cell +
				geometry->cell( particles->geo_state(pidx) );
	atomicAdd( &tally[tally_idx], particles->wt(pidx) );
    }
}

//---------------------------------------------------------------------------//
// Finalize the tally.
__global__ void finalize_kernel( const int num_batch,
				 const int num_cell,
				 const int total_num_particle,
				 double* tally )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < num_batch * num_cell ) 
    {
	tally[idx] *= num_batch / total_num_particle;
    }
}

//---------------------------------------------------------------------------//
// Calculate the first and second moments of the tally.
__global__ void moments_kernel( const int num_batch,
				const int num_cell,
				const double* tally,
				double* first_moment,
				double* second_moment )
{
    std::size_t cell_idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( cell_id < num_cell )
    {
	int tally_idx = 0;

	// Calculate the first moment.
	first_moment[cell_idx] = 0.0;
	for ( int b = 0; b < num_batch; ++b )
	{
	    tally_idx = b * num_cells + cell_idx;
	    first_moment[cell_idx] += tally[ tally_idx ];
	}
	first_moment[cell_idx] /= num_batch;

	// Calculate the second moment.
	second_moment[cells_idx] = 0.0;
	for ( int b = 0; b < num_batch; ++b )
	{
	    tally_idx = b * num_cells + cell_idx;
	    second_moment[cell_idx] += 
		tally[ tally_idx ] * tally[ tally_idx ]
		- first_moment[cell_idx] * first_moment[cell_idx];
	}
	second_moment[cell_idx] /= num_batch * (num_batch - 1);
    }
}

//---------------------------------------------------------------------------//
// HOST API
//---------------------------------------------------------------------------//
// Constructor.
template <class Geometry>
Collision_Tally<Geometry>::Collision_Tally( 
    const cuda::Shared_Device_Ptr<Geometry>& geometry, 
    const int num_batch )
    : d_geometry( d_geometry )
    , d_num_batch( num_batch )
    , d_num_cells( d_geometry.get_host_ptr()->num_cells() )
{
    // Allocate the tally.
    std::size_t size = d_num_batch * d_num_cells;
    std::size_t mem_size = size * sizeof(double);
    CudaCall( cudaMalloc( (void**) &d_tally, mem_size ) );
    
    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    // Initialize the tally to zero.
    init_tally_kernel<<<num_blocks,threads_per_block>>>( d_size, d_tally );
}
    
//---------------------------------------------------------------------------//
// Destructor.
template <class Geometry>
Collision_Tally<Geometry>::~Collision_Tally()
{
    cudaFree( d_tally );
}

//---------------------------------------------------------------------------//
// Tally the particles in a vector.
template <class Geometry>
void Collision_Tally<Geometry>::tally( 
    const cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles )
{
    // Get the particles that have had a collision.
    std::size_t start_idx = 0;
    std::size_t num_particle = 0;
    particles.get_host_ptr()->get_event_particles( profugus::Events::COLLISION,
						   start_index,
						   num_particle );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = num_particle / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Tally the particles.
    tally_kernel<<<num_blocks,threads_per_block>>>( d_geometry.get_device_ptr(),
						    particles.get_device_ptr(),
						    start_idx,
						    num_particle,
						    d_num_batch,
						    d_num_cells,
						    d_tally );
}

//---------------------------------------------------------------------------//
// Finalize the tally.
template <class Geometry>
void Collision_Tally<Geometry>::finalize( const std::size_t total_num_particle )
{
    // Get CUDA launch parameters.
    int size = d_num_batch * d_num_cells;
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = size / threads_per_block;
    if ( size % threads_per_block > 0 ) ++num_blocks;

    // Finalize the tally.
    finalize_kernel<<<num_blocks,threads_per_block>>>( d_num_batch,
						       d_num_cell,
						       total_num_particle,
						       d_tally );
}

//---------------------------------------------------------------------------//
// Copy the tally moments to the host.
template <class Geometry>
void Collision_Tally<Geometry>::copy_moments_to_host( 
    Teuchos::Array<double>& first_moment,
    Teuchos::Array<double>& second_moment )
{
    // Allocate moments on device.
    std::size_t mem_size = d_num_cells * sizeof(double);
    double* device_first_moment;
    CudaCall( cudaMalloc( (void**) &device_first_moment, mem_size ) );
    double* device_second_moment;
    CudaCall( cudaMalloc( (void**) &device_second_moment, mem_size ) );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = size / threads_per_block;
    if ( d_num_cells % threads_per_block > 0 ) ++num_blocks;

    // Calculate moments on device.
    moments_kernel<<<num_blocks,threads_per_block>>>( d_num_batch,
						      d_num_cell,
						      d_tally,
						      device_first_moment,
						      device_second_moment );

    // Copy moments to host.
    first_moment.resize( d_num_cells );
    CudaCall( cudaMemcpy( first_moment.getRawPtr(), device_first_moment,
			  mem_size, cudaMemcpyDeviceToHost) );
    second_moment.resize( d_num_cells );
    CudaCall( cudaMemcpy( second_moment.getRawPtr(), device_second_moment,
			  mem_size, cudaMemcpyDeviceToHost) );

    // Free device moments.
    cudaFree( device_first_moment );
    cudaFree( device_second_moment );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // end cuda_mc_Collision_Tally_t_cuh

//---------------------------------------------------------------------------//
//                 end of Collision_Tally.t.cuh
//---------------------------------------------------------------------------//
