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

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Cell_Tally<Geometry>::Cell_Tally(SDP_Geometry geometry,
                                 SDP_Physics  physics)
    : d_geometry_shared(geometry)
{
    REQUIRE(d_geometry_shared.get_host_ptr());
    REQUIRE(d_geometry_shared.get_device_ptr());
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

    d_host_cells = cells;

    d_num_cells = d_host_cells.size();
    REQUIRE( d_num_cells > 0 );

    d_cell_vec.resize(d_num_cells);
    d_tally_vec.resize(d_num_cells,0.0);

    // Extract pointers for on-device access
    d_cells  = d_cell_vec.data().get();
    d_tally  = d_tally_vec.data().get();

    // Copy cells to device
    cudaMemcpy( d_cells, &d_host_cells[0], d_num_cells*sizeof(int),
        cudaMemcpyHostToDevice);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Compute volumes for tally cells
    d_host_volumes.resize(d_num_cells);
    const auto &geom_host = d_geometry_shared.get_host_ptr();
    const auto &all_volumes = geom_host->volumes();
    for( int ind = 0; ind < d_num_cells; ++ind )
    {
        int cell = d_host_cells[ind];
        d_host_volumes[ind] = all_volumes[cell];
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing of tally
 */
template <class Geometry>
void Cell_Tally<Geometry>::finalize(double                num_particles,
                                    Cell_Tally<Geometry> *cell_tally_dev)
{
    REQUIRE(num_particles > 1);

    // Copy results to host to normalize tally results
    d_host_tally.resize(d_num_cells);
    cudaMemcpy( &d_host_tally[0], d_tally, d_num_cells * sizeof(double),
                cudaMemcpyDeviceToHost );

    // Global reduction on tally results
    profugus::global_sum(&d_host_tally[0],d_num_cells);

    // Now normalize tally
    for( int ind = 0; ind < d_num_cells; ++ind )
    {
        double norm_factor = 1.0 / (d_host_volumes[ind] * num_particles);
        d_host_tally[ind] *= norm_factor;
    }

    std::cout << "Cell tally: ";
    for( const auto &t : d_host_tally )
        std::cout << t << " ";
    std::cout << std::endl;
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
