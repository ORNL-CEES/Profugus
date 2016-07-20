//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Cell_Tally.t.cu
 * \author Stuart Slattery
 * \brief  Cell class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Cell_Tally_t_cuh
#define cuda_mc_Cell_Tally_t_cuh

#include "cuda_utils/Hardware.hh"
#include "cuda_utils/Memory.cuh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Atomic_Add.cuh"

#include "mc/Definitions.hh"

#include "Cell_Tally.hh"

#include <cuda_runtime.h>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Initialize the tally to zero
__global__ void init_tally_kernel( const int size, double* tally )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) tally[idx] = 0.0;
}

//---------------------------------------------------------------------------//
// Tally particles.
template<class Geometry>
__global__ void tally_kernel( const Geometry* geometry,
			      const Particle_Vector<Geometry>* particles,
			      const int num_collision,
			      const int num_boundary,
			      const int num_batch,
			      const int num_cell,
			      double* tally )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Get the event indices.
    int* collision_indices = particles->event_indices( events::COLLISION );
    int* boundary_indices = particles->event_indices( events::BOUNDARY );

    if ( idx < num_collision + num_boundary )
    {
	// Get the particle index.
	int pidx = ( idx < num_collision )
		   ? collision_indices[idx]
		   : boundary_indices[idx - num_collision];
    
	// Accumulate the particle in its batch and cell.
	REQUIRE( particles->alive(pidx) );
	int tally_idx = particles->batch( pidx ) * num_cell +
				geometry->cell( particles->geo_state(pidx) );
	CHECK( tally_idx < num_batch * num_cell );
	cuda::Atomic_Add<cuda::arch::Device>::fetch_add( 
	    &tally[tally_idx], particles->wt(pidx) * particles->step(pidx) );
    }
}

//---------------------------------------------------------------------------//
// Finalize the tally.
__global__ void finalize_kernel( const int num_batch,
				 const int num_cell,
				 const double num_particles,
				 double* tally )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < num_batch * num_cell ) 
    {
	tally[idx] = (tally[idx] * num_batch) / num_particles;
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
    int cell_idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( cell_idx < num_cell )
    {
	int tally_idx = 0;

	// Calculate the first moment.
	first_moment[cell_idx] = 0.0;
	for ( int b = 0; b < num_batch; ++b )
	{
	    tally_idx = b * num_cell + cell_idx;
	    CHECK( tally_idx < num_batch * num_cell );
	    first_moment[cell_idx] += tally[ tally_idx ];
	}
	first_moment[cell_idx] /= num_batch;

	// Calculate the second moment.
	second_moment[cell_idx] = 0.0;
	for ( int b = 0; b < num_batch; ++b )
	{
	    tally_idx = b * num_cell + cell_idx;
	    CHECK( tally_idx < num_batch * num_cell );
	    second_moment[cell_idx] += 
		tally[ tally_idx ] * tally[ tally_idx ] -
		first_moment[cell_idx] * first_moment[cell_idx];
	}
	second_moment[cell_idx] /= num_batch * (num_batch - 1);
    }
}

//---------------------------------------------------------------------------//
// HOST API
//---------------------------------------------------------------------------//
// Constructor.
template <class Geometry>
Cell_Tally<Geometry>::Cell_Tally( 
    const cuda::Shared_Device_Ptr<Geometry>& geometry, 
    const int num_batch )
    : d_geometry( geometry )
    , d_num_batch( num_batch )
    , d_num_cells( d_geometry.get_host_ptr()->num_cells() )
{
    // Allocate the tally.
    int size = d_num_batch * d_num_cells;
    cuda::memory::Malloc( d_tally, size );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = size / threads_per_block;
    if ( size % threads_per_block > 0 ) ++num_blocks;

    // Initialize the tally to zero.
    init_tally_kernel<<<num_blocks,threads_per_block>>>( size, d_tally );
}
    
//---------------------------------------------------------------------------//
// Destructor.
template <class Geometry>
Cell_Tally<Geometry>::~Cell_Tally()
{
    cuda::memory::Free( d_tally );
}

//---------------------------------------------------------------------------//
// Tally the particles in a vector.
template <class Geometry>
void Cell_Tally<Geometry>::accumulate( 
    const cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles )
{
    // Get the particles that just had a collision.
    int num_collision =
        particles.get_host_ptr()->get_event_size( events::COLLISION );

    // Get the particles that just hit a boundary.
    int num_boundary = 
        particles.get_host_ptr()->get_event_size( events::BOUNDARY );

    // Calculate the launch parameters.
    int num_particle = num_collision + num_boundary;
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = num_particle / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Tally the particles.
    tally_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>( 
        d_geometry.get_device_ptr(),
        particles.get_device_ptr(),
        num_collision,
        num_boundary,
        d_num_batch,
        d_num_cells,
        d_tally );

    // Synchronize after tally.
    d_stream.synchronize();
}

//---------------------------------------------------------------------------//
// Finalize the tally.
template <class Geometry>
void Cell_Tally<Geometry>::finalize( double num_particles )
{
    // Get CUDA launch parameters.
    int size = d_num_batch * d_num_cells;
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = size / threads_per_block;
    if ( size % threads_per_block > 0 ) ++num_blocks;

    // Finalize the tally.
    finalize_kernel<<<num_blocks,threads_per_block>>>( d_num_batch,
						       d_num_cells,
						       num_particles,
						       d_tally );
}

//---------------------------------------------------------------------------//
// Copy the tally moments to the host.
template <class Geometry>
void Cell_Tally<Geometry>::copy_moments_to_host( 
    Teuchos::Array<double>& first_moment,
    Teuchos::Array<double>& second_moment ) const
{
    // Allocate moments on device.
    double* device_first_moment = NULL;
    cuda::memory::Malloc( device_first_moment, d_num_cells );
    double* device_second_moment = NULL;
    cuda::memory::Malloc( device_second_moment, d_num_cells );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::default_block_size();
    unsigned int num_blocks = d_num_cells / threads_per_block;
    if ( d_num_cells % threads_per_block > 0 ) ++num_blocks;

    // Calculate moments on device.
    moments_kernel<<<num_blocks,threads_per_block>>>( d_num_batch,
						      d_num_cells,
						      d_tally,
						      device_first_moment,
						      device_second_moment );

    // Copy moments to host.
    first_moment.resize( d_num_cells );
    cuda::memory::Copy_To_Host( 
	first_moment.getRawPtr(), device_first_moment, d_num_cells );
    second_moment.resize( d_num_cells );
    cuda::memory::Copy_To_Host( 
	second_moment.getRawPtr(), device_second_moment, d_num_cells );

    // Free device moments.
    cuda::memory::Free( device_first_moment );
    cuda::memory::Free( device_second_moment );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // end cuda_mc_Cell_Tally_t_cuh

//---------------------------------------------------------------------------//
//                 end of Cell_Tally.t.cuh
//---------------------------------------------------------------------------//
