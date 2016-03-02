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

// Functor to compute cell volumes on device
template <class Geometry>
class Volume_Functor
{
  public:

      typedef cuda::arch::Device                 Arch_t;
      typedef cuda::Device_Vector<Arch_t,int>    Dev_Int;
      typedef cuda::Device_Vector<Arch_t,double> Dev_Dbl;

      Volume_Functor( Geometry *geom, Dev_Int &cells, Dev_Dbl &volumes )
          : d_geometry(geom)
          , d_cells(cells)
          , d_volumes(volumes)
      {
          ENSURE( d_cells.size() == d_volumes.size() );
      }

      // Compute volume with indirection from cell list
      __device__ void operator()(std::size_t ind)
      {
          int cell = d_cells[ind];
          d_volumes[ind] = d_geometry->volume(cell);
      }

  private:

      Geometry *d_geometry;
      Dev_Int   d_cells;
      Dev_Dbl   d_volumes;
};

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Cell_Tally<Geometry>::Cell_Tally(SDP_Geometry geometry,
                                 SDP_Physics  physics)
    : d_tally(0)
    , d_cells(0)
    , d_volumes(0)
{
    REQUIRE(geometry.get_host_ptr());
    REQUIRE(geometry.get_device_ptr());
    REQUIRE(physics.get_host_ptr());
    REQUIRE(physics.get_device_ptr());

    d_geometry = geometry.get_device_ptr();
    d_physics  = physics.get_device_ptr();
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Add a list of cells to tally.
 */
template <class Geometry>
void Cell_Tally<Geometry>::set_cells(const std::vector<int> &cells)
{
    REQUIRE( std::is_sorted(cells.begin(),cells.end()) );

    auto num_cells = cells.size();

    // Make local device vectors for cells and tally result,
    // then swap with members
    Dev_Int dev_cells(profugus::make_view(cells));
    cuda::swap(d_cells,dev_cells);

    std::vector<double> zeros(num_cells,0.0);
    Dev_Dbl dev_data(profugus::make_view(zeros));
    cuda::swap(d_tally,dev_data);

    // Compute volumes for tally cells
    Dev_Dbl dev_volumes(num_cells);
    cuda::Launch_Args<Arch_t> launch_args;
    launch_args.set_num_elements(num_cells);

    Volume_Functor<Geometry> volume_func( d_geometry, d_cells, dev_volumes );
    cuda::parallel_launch( volume_func, launch_args );

    cuda::swap( d_volumes, dev_volumes );
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing on first and second moments.
 */
template <class Geometry>
void Cell_Tally<Geometry>::finalize(double num_particles)
{
    REQUIRE(num_particles > 1);

    int num_cells = d_tally.size();
    REQUIRE( d_cells.size()   == num_cells );
    REQUIRE( d_volumes.size() == num_cells );

    // Copy results to host to normalize tally results
    cuda::Host_Vector<double> host_tally(num_cells);
    d_tally.to_host( host_tally );

    // Global reduction on tally results
    profugus::global_sum(&host_tally[0],host_tally.size());

    // Get host copy of cells
    cuda::Host_Vector<double> host_volumes(num_cells);
    d_volumes.to_host( host_volumes );

    // Now normalize tally
    for( int ind = 0; ind < num_cells; ++ind )
    {
        double norm_factor = 1.0 / (host_volumes[ind] * num_particles);
        host_tally[ind] *= norm_factor;
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Clear/re-initialize all tally values between solves.
 */
template <class Geometry>
void Cell_Tally<Geometry>::reset()
{
    // Clear all current tally results (keep cell list)
    std::vector<double> z(d_tally.size(),0.0);
    d_tally.assign( profugus::make_view(z) );
}


} // end namespace cuda_mc

#endif // MC_cuda_mc_Cell_Tally_t_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Cell_Tally.t.cuh
//---------------------------------------------------------------------------//
