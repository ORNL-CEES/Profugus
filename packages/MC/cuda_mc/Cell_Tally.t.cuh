//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Cell_Tally.t.cuh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Cell_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Cell_Tally_t_cuh
#define MC_cuda_mc_Cell_Tally_t_cuh

#include <cmath>
#include <algorithm>

#include "Utils/comm/global.hh"
#include "Utils/utils/View_Field.hh"

#include "cuda_utils/Launch_Args.t.cuh"
#include "cuda_utils/Host_Vector.hh"

#include "Cell_Tally.cuh"

namespace cuda_mc
{

// Function to compute cell volumes on device
template <class Geometry>
__global__ void compute_volumes( Geometry  *geom,
                                 const int *cells,
                                 double    *volumes,
                                 int        num_cells )
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_cells )
        volumes[tid] = geom->volume(cells[tid]);
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Cell_Tally<Geometry>::Cell_Tally(SDP_Geometry geometry,
                                 SDP_Physics  physics)
{
    REQUIRE(geometry.get_host_ptr());
    REQUIRE(geometry.get_device_ptr());
    REQUIRE(physics.get_host_ptr());
    REQUIRE(physics.get_device_ptr());

    d_geometry = geometry.get_device_ptr();
    d_physics  = physics.get_device_ptr();
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing of tally
 */
template <class Geometry>
void Cell_Tally<Geometry>::set_cells(const std::vector<int> &cells)
{
    REQUIRE( std::is_sorted(cells.begin(),cells.end()) );

    d_num_cells = cells.size();
    REQUIRE( d_num_cells > 0 );

    cudaMalloc( (void**)&d_cells, d_num_cells * sizeof(int) );
    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaMalloc( (void**)&d_volumes, d_num_cells * sizeof(double) );
    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaMalloc( (void**)&d_tally, d_num_cells * sizeof(double) );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy cells to device
    cudaMemcpy( d_cells, &cells[0], d_num_cells*sizeof(int),
        cudaMemcpyHostToDevice);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Initialize tally to zero
    std::vector<double> zeros(d_num_cells,0.0);
    cudaMemcpy( d_tally, &zeros[0], d_num_cells*sizeof(double),
        cudaMemcpyHostToDevice);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Compute volumes for tally cells
    int block_size = std::min(256,d_num_cells);
    compute_volumes<<<1,block_size>>>(d_geometry, d_cells,
                                      d_volumes,  d_num_cells);
    REQUIRE( cudaGetLastError() == cudaSuccess );

}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing of tally
 */
template <class Geometry>
void Cell_Tally<Geometry>::finalize(double num_particles)
{
    REQUIRE(num_particles > 1);

    // Copy results to host to normalize tally results
    std::vector<double> host_tally(d_num_cells);
    cudaMemcpy( &host_tally[0], d_tally, d_num_cells*sizeof(double),
                cudaMemcpyDeviceToHost );

    // Global reduction on tally results
    profugus::global_sum(&host_tally[0],host_tally.size());

    // Get host copy of cells
    std::vector<double> host_volumes(d_num_cells);
    cudaMemcpy( &host_volumes[0], d_volumes, d_num_cells*sizeof(double),
                cudaMemcpyDeviceToHost );

    // Now normalize tally
    for( int ind = 0; ind < d_num_cells; ++ind )
    {
        double norm_factor = 1.0 / (host_volumes[ind] * num_particles);
        host_tally[ind] *= norm_factor;
    }

    // Copy results back to device
    cudaMemcpy( d_tally, &host_tally[0], d_num_cells*sizeof(double),
                cudaMemcpyHostToDevice);
}

//---------------------------------------------------------------------------//
/*
 * \brief Clear/re-initialize all tally values between solves.
 */
template <class Geometry>
void Cell_Tally<Geometry>::reset()
{
    // Clear all current tally results (keep cell list)
    std::vector<double> z(d_num_cells,0.0);
    cudaMemcpy( d_tally, &z[0], d_num_cells*sizeof(double),
                cudaMemcpyHostToDevice );
}


} // end namespace cuda_mc

#endif // MC_cuda_mc_Cell_Tally_t_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Cell_Tally.t.cuh
//---------------------------------------------------------------------------//
