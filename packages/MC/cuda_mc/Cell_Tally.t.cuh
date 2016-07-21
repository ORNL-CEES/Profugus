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
#include "Utils/utils/Serial_HDF5_Writer.hh"

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

    // Initilialize batch count
    d_tally_sum.resize(d_num_cells,0.0);
    d_tally_sum_sq.resize(d_num_cells,0.0);
}

//---------------------------------------------------------------------------//
/*
 * \brief End batch for statistics
 */
template <class Geometry>
void Cell_Tally<Geometry>::end_batch(double num_particles)
{
    // Copy results to host to normalize tally results
    d_host_tally.resize(d_num_cells);
    cudaMemcpy( &d_host_tally[0], d_tally, d_num_cells * sizeof(double),
                cudaMemcpyDeviceToHost );

    // Global reduction on tally results
    profugus::global_sum(&d_host_tally[0],d_num_cells);

    for (int cell = 0; cell < d_num_cells; ++cell)
    {
        auto val = d_host_tally[cell];
        d_tally_sum[cell]    += val;
        d_tally_sum_sq[cell] += val*val / num_particles;
    }

    d_batch_np.push_back(num_particles);

    // Zero out tally for next cycle
    thrust::fill(d_host_tally.begin(),d_host_tally.end(),0.0);
    d_tally_vec = d_host_tally;
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing of tally
 */
template <class Geometry>
void Cell_Tally<Geometry>::finalize(double num_particles)
{
    REQUIRE(num_particles > 0);

    // If no batches have been processed, make sure end_batch is called
    if (d_batch_np.empty())
        end_batch(num_particles);

    // Make sure number of particles at finalization matches sum of batches
    double sum_np = std::accumulate(d_batch_np.begin(),d_batch_np.end(),0.0);
    REQUIRE(profugus::soft_equiv(num_particles,sum_np));

    double num_batches = static_cast<double>(d_batch_np.size());

    d_host_std_dev.resize(d_num_cells,0.0);

    // Now compute std dev and normalize tally
    for (int cell = 0; cell < d_num_cells; ++cell)
    {
        double tally_mean = d_tally_sum[cell] / num_particles;
        double var = (d_tally_sum_sq[cell] / num_particles) -
                     tally_mean * tally_mean;

        if (var > 0.0)
            var /= num_batches;

        // Normalize by volume
        double vol = d_host_volumes[cell];
        d_host_tally[cell] = tally_mean / vol;
        d_host_std_dev[cell] = std::sqrt(var) / vol;
    }

#ifdef USE_HDF5
    // Open file for writing
    std::string filename = "cell_tally.h5";
    profugus::Serial_HDF5_Writer writer;
    writer.open(filename);

    // Write data
    writer.write("cells",d_host_cells);
    writer.write("flux_mean",d_host_tally);
    writer.write("flux_std_dev",d_host_std_dev);

    writer.close();
#endif
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

    std::fill(d_tally_sum.begin(),d_tally_sum.end(),0.0);
    std::fill(d_tally_sum_sq.begin(),d_tally_sum_sq.end(),0.0);
    d_batch_np.clear();
}


} // end namespace cuda_mc

#endif // MC_cuda_mc_Cell_Tally_t_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Cell_Tally.t.cuh
//---------------------------------------------------------------------------//
