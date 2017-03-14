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

#include "utils/Serial_HDF5_Writer.hh"

#include "comm/global.hh"

#include "mc/Definitions.hh"

#include "Cell_Tally.hh"

#include <cuda_runtime.h>

#include <cmath>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//

namespace
{

// Initialize the tally to zero
__global__ void init_tally_kernel( const int size, double* tally )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < size ) tally[idx] = 0.0;
}

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
    int collision_start = particles->event_lower_bound( events::COLLISION );
    int boundary_start = particles->event_lower_bound( events::BOUNDARY );

    if ( idx < num_collision + num_boundary )
    {
	// Get the particle index.
	int pidx = ( idx < num_collision )
		   ? idx + collision_start
		   : idx - num_collision + boundary_start;
    
	// Accumulate the particle in its batch and cell.
	REQUIRE( particles->alive(pidx) );
	int tally_idx = particles->batch( pidx ) * num_cell +
				geometry->cell( particles->geo_state(pidx) );
	CHECK( tally_idx < num_batch * num_cell );
	cuda_utils::Atomic_Add<cuda_utils::arch::Device>::fetch_add( 
	    &tally[tally_idx], particles->wt(pidx) * particles->step(pidx) );
    }
}

//---------------------------------------------------------------------------//
// HOST API
//---------------------------------------------------------------------------//
// Constructor.
template <class Geometry>
Cell_Tally<Geometry>::Cell_Tally( 
    RCP_Std_DB db,
    const cuda_utils::Shared_Device_Ptr<Geometry>& geometry, 
    const int num_batch )
    : d_geometry( geometry )
    , d_num_batch( num_batch )
    , d_num_cells( d_geometry.get_host_ptr()->num_cells() )
    , d_db(db)
{
    // Allocate the tally.
    int size = d_num_batch * d_num_cells;
    cuda_utils::memory::Malloc( d_tally, size );

    // Get CUDA launch parameters.
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
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
    cuda_utils::memory::Free( d_tally );
}

//---------------------------------------------------------------------------//
// Tally the particles in a vector.
template <class Geometry>
void Cell_Tally<Geometry>::accumulate( 
    const cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles )
{
    // Get the particles that just had a collision.
    int num_collision =
        particles.get_host_ptr()->get_event_size( events::COLLISION );

    // Get the particles that just hit a boundary.
    int num_boundary = 
        particles.get_host_ptr()->get_event_size( events::BOUNDARY );

    // Calculate the launch parameters.
    int num_particle = num_collision + num_boundary;
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
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
    // Copy the local tally to the host.
    int size = d_num_batch * d_num_cells;
    std::vector<double> batch_results( size );
    cuda_utils::memory::Copy_To_Host( batch_results.data(), d_tally, size );

    // Calculate the batch results.
    int total_num_batch = (profugus::nodes() * d_num_batch);
    int particles_per_batch = num_particles / total_num_batch;
    for ( int c = 0; c < d_num_cells; ++c )
    {
        for ( int b = 0; b < d_num_batch; ++b )
        {
            batch_results[ b * d_num_cells + c ] /= 
              particles_per_batch * 
              d_geometry.get_host_ptr()->mesh().volume_host( c );
        }
    }    

    // Calculate the first moment.
    d_first.resize( d_num_cells, 0.0 );
    for ( int c = 0; c < d_num_cells; ++c )
    {
        for ( int b = 0; b < d_num_batch; ++b )
        {
            d_first[ c ] += batch_results[ b * d_num_cells + c ];
        }
    }
    profugus::global_sum(d_first.data(), d_first.size());
    for ( int c = 0; c < d_num_cells; ++c )
      {
        d_first[ c ] /= total_num_batch;
      }      

    // Calculate the second moment (standard deviation).
    d_second.resize( d_num_cells, 0.0 );
    for ( int c = 0; c < d_num_cells; ++c )
    {
        for ( int b = 0; b < d_num_batch; ++b )
        {
            d_second[ c ] += batch_results[ b * d_num_cells + c ] *
                           batch_results[ b * d_num_cells + c ];
        }
    }
    profugus::global_sum(d_second.data(), d_second.size());
    for ( int c = 0; c < d_num_cells; ++c )
    {
      d_second[ c ] /= total_num_batch;
      d_second[ c ] -= d_first[c]*d_first[c];
      d_second[ c ] /= total_num_batch - 1;
      d_second[ c ] = std::sqrt( d_second[ c ] );
    }
    
    
    // Determine permutation vector for sorted cells
    std::vector<int>    cells(d_num_cells,  0);
    std::vector<std::pair<int,int>> sort_vec;
    for (int i = 0; i < cells.size(); ++i)
      {
        cells[i] = i;
        sort_vec.push_back({cells[i],i});
      }

    std::sort(sort_vec.begin(), sort_vec.end(),
        [](const std::pair<int,int> &lhs, const std::pair<int,int> &rhs)
        { return lhs.first < rhs.first; } );

    std::vector<int> cell_map(cells.size());
    for (int i = 0; i < sort_vec.size(); ++i)
        cell_map[i] = sort_vec[i].second;

    // Reorder vectors
    {
        std::vector<int>    tmp_cells  = cells;
        std::vector<double> tmp_first  = d_first;
        std::vector<double> tmp_second = d_second;
        for (int i = 0; i < cell_map.size(); ++i)
        {
            int ind = cell_map[i];
            cells[i]  = tmp_cells[ind];
            d_first[i]  = tmp_first[ind];
            d_second[i] = tmp_second[ind];
        }
    }

#ifdef USE_HDF5
    REQUIRE( d_db->isType<std::string>("problem_name") );
    std::string filename = d_db->get<std::string>("problem_name") + "_flux.h5";

    profugus::Serial_HDF5_Writer writer;
    writer.open(filename);
    writer.write("cells",cells);
    writer.write("flux_mean",d_first);
    writer.write("flux_std_dev",d_second);
    writer.close();
#endif
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // end cuda_mc_Cell_Tally_t_cuh

//---------------------------------------------------------------------------//
//                 end of Cell_Tally.t.cuh
//---------------------------------------------------------------------------//
