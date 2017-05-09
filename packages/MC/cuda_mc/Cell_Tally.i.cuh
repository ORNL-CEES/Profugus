//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Cell_Tally.i.cuh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Cell_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Cell_Tally_i_cuh
#define MC_cuda_mc_Cell_Tally_i_cuh

#include "cuda_utils/Utility_Functions.hh"
#include "Cell_Tally.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*
 * \brief Track particle and tally..
 */
template <class Geometry>
__device__ void Cell_Tally<Geometry>::accumulate(
        double                   step,
        int                      pid,
        const Particle_Vector_t *particles)
{
    DEVICE_REQUIRE( step >= 0.0 );
    DEVICE_REQUIRE( particles->alive(pid) );

    // Get the cell index
    int cell = d_geometry->cell(particles->geo_states(),pid);

    int ind = cuda::utility::lower_bound( d_cells, 
                                          d_cells+d_num_cells,
                                          cell ) - d_cells;

    // See if valid index was returned and if it is an exact match
    // for the particle's cell
    if( (ind < d_num_cells) && (d_cells[ind] == cell))
    {
        cuda::utility::atomic_add_double(&d_tally[ind],
                                         particles->wt(pid) * step );
    }
}

} // namespace cuda_mc

#endif // MC_cuda_mc_Cell_Tally_i_cuh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.i.hh
//---------------------------------------------------------------------------//
