//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Geometry.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 28 16:37:36 2014
 * \brief  Geometry member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "core/geometry/Definitions.hh"
#include "Geometry.hh"

namespace gpu
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Geometry::Geometry(int    N,
                   double d)
    : d_edges(3, std::vector<double>(N+1, 0.0))
{
    // build the grid edges
    for (int n = 1; n < N+1; ++n)
    {
        d_edges[0][n] = d_edges[0][n-1] + d;
        d_edges[1][n] = d_edges[1][n-1] + d;
        d_edges[2][n] = d_edges[2][n-1] + d;
    }

    d_x = &d_edges[0][0];
    d_y = &d_edges[1][0];
    d_z = &d_edges[2][0];

    std::fill(std::begin(d_N), std::end(d_N), N);

#pragma acc enter data copyin(this)
#pragma acc enter data copyin(d_x[0:N+1], d_y[0:N+1], d_z[0:N+1])
#pragma acc enter data copyin(d_N[0:3])
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Geometry::~Geometry()
{
#pragma acc exit data delete(this)
#pragma acc exit data delete(d_x, d_y, d_z, d_N)
}

//---------------------------------------------------------------------------//
// GEOMETRY FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate distance to the next cell.
 */
double Geometry::distance_to_boundary(Geo_State_t& state) const
{
    using def::I; using def::J; using def::K;

    // large double for comparisons
    double max_d = 1.0e24;

    // min distance to boundary; initializing test_dist to a large number before
    // each surface check implicitly handles the case where omega[dir] == 0.0
    double test_dist = max_d;
    int    test_next = state.ijk[I];

    // initialize the next surface
    state.next_ijk = state.ijk;

    // unrolled check

    // X SURFACE
    if (state.dir[I] > 0.0 && state.ijk[I] < d_N[I])
    {
        test_dist = (d_x[state.ijk[I]+1] - state.pos[I]) / state.dir[I];
        test_next = state.ijk[I] + 1;
    }
    else if (state.dir[I] < 0.0 && state.ijk[I] > -1)
    {
        test_dist = (d_x[state.ijk[I]] - state.pos[I]) / state.dir[I];
        test_next = state.ijk[I] - 1;
    }

    // initialize to x distance
    state.next_dist   = test_dist;
    state.next_ijk[I] = test_next;

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = max_d;

    // Y SURFACE
    if (state.dir[J] > 0.0 && state.ijk[J] < d_N[J])
    {
        test_dist = (d_y[state.ijk[J]+1] - state.pos[J]) / state.dir[J];
        test_next = state.ijk[J] + 1;
    }
    else if (state.dir[J] < 0.0 && state.ijk[J] > -1)
    {
        test_dist = (d_y[state.ijk[J]] - state.pos[J]) / state.dir[J];
        test_next = state.ijk[J] - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist   = test_dist;
        state.next_ijk[I] = state.ijk[I];
        state.next_ijk[J] = test_next;
    }

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = max_d;

    // Z SURFACE
    if (state.dir[K] > 0.0 && state.ijk[K] < d_N[K])
    {
        test_dist = (d_z[state.ijk[K]+1] - state.pos[K]) / state.dir[K];
        test_next = state.ijk[K] + 1;
    }
    else if (state.dir[K] < 0.0 && state.ijk[K] > -1)
    {
        test_dist = (d_z[state.ijk[K]] - state.pos[K]) / state.dir[K];
        test_next = state.ijk[K] - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist   = test_dist;
        state.next_ijk[I] = state.ijk[I];
        state.next_ijk[J] = state.ijk[J];
        state.next_ijk[K] = test_next;
    }

    return state.next_dist;
}

} // end namespace gpu

//---------------------------------------------------------------------------//
//                 end of Geometry.cc
//---------------------------------------------------------------------------//
