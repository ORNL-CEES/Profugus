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

#include "cuda_utility/Utility_Functions.hh"
#include "Cell_Tally.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*
 * \brief Track particle and tally..
 */
template <class Geometry>
void Cell_Tally<Geometry>::accumulate(double            step,
                                      const Particle_t &p)
{
    REQUIRE( step >= 0.0 );
    REQUIRE( p.alive() );

    // Get the cell index
    int cell = d_geometry->cell(p.geo_state());

    int ind = cuda::utility::lower_bound( &d_cells[0], 
                                          &d_cells[0]+d_cells.size(),
                                           cell ) - &d_cells[0];

    CHECK( ind < d_tally.size() );

    cuda::atomic_add_double( &d_tally[ind], p.wt() * step );
}


#endif // MC_cuda_mc_Cell_Tally_i_cuh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.i.hh
//---------------------------------------------------------------------------//
