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

#include "Geometry.hh"
#include "Change_Direction.hh"

namespace acc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Geometry::Geometry(int                     N,
                   double                  d,
                   const std::vector<int> &matids,
                   const int              *bnds)
    : d_edges(3, std::vector<double>(N+1, 0.0))
    , d_matids(matids)
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
    d_m = &d_matids[0];

    int nc = num_cells();

    std::fill(std::begin(d_N), std::end(d_N), N);
    std::copy(bnds, bnds+6, std::begin(d_b));
    std::copy(bnds, bnds+6, std::begin(d_bnds));

#pragma acc enter data pcopyin(this)
#pragma acc enter data pcopyin(d_x[0:N+1], d_y[0:N+1], d_z[0:N+1])
#pragma acc enter data pcopyin(d_N[0:3])
#pragma acc enter data pcopyin(d_b[0:6], d_m[0:nc])
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Geometry::~Geometry()
{
#pragma acc exit data delete(d_x, d_y, d_z, d_N, d_b, d_m)
#pragma acc exit data delete(this)
}

//---------------------------------------------------------------------------//
// GEOMETRY FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize particle
 *
 * This finds the closest ijk coords on each axis to the location. A particle
 * can be born "outside" and have ijk extents that are outside [0,N) .
 */
void Geometry::initialize(const double   *r,
                          const double   *direction,
                          Geometry_State &state)
{

    // Set struct attributes
    state.pos[0] = r[0];
    state.pos[1] = r[1];
    state.pos[2] = r[2];

    state.dir[0] = direction[0];
    state.dir[1] = direction[1];
    state.dir[2] = direction[2];

    state.ijk[0] = state.pos[0] / (d_x[1] - d_x[0]);
    state.ijk[1] = state.pos[1] / (d_y[1] - d_y[0]);
    state.ijk[2] = state.pos[2] / (d_z[1] - d_z[0]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate distance to the next cell.
 */
double Geometry::distance_to_boundary(Geometry_State& state) const
{
    using def::I; using def::J; using def::K;

    // large double for comparisons
    double max_d = 1.0e24;

    // min distance to boundary; initializing test_dist to a large number before
    // each surface check implicitly handles the case where omega[dir] == 0.0
    double test_dist = max_d;
    int    test_next = state.ijk[I];

    // initialize the next surface
    state.next_ijk[0] = state.ijk[0];
    state.next_ijk[1] = state.ijk[1];
    state.next_ijk[2] = state.ijk[2];

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

//---------------------------------------------------------------------------//
/*!
 * \brief Change the direction through an angle.
 */
void Geometry::change_direction(
        double          costheta,
        double          phi,
        Geometry_State &state) const
{
    ::acc::change_direction(costheta, phi, state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the boundary state.
 */
int Geometry::boundary_state(const Geometry_State &state) const
{
    using def::I; using def::J; using def::K;

    if (d_b[0] && state.next_ijk[0] == -1)
    {
        return 2;
    }
    else if (d_b[1] && state.next_ijk[0] == d_N[0])
    {
        return 2;
    }

    if (d_b[2] && state.next_ijk[1] == -1)
    {
        return 2;
    }
    else if (d_b[3] && state.next_ijk[1] == d_N[1])
    {
        return 2;
    }

    if (d_b[4] && state.next_ijk[2] == -1)
    {
        return 2;
    }
    else if (d_b[5] && state.next_ijk[2] == d_N[2])
    {
        return 2;
    }

    if ((state.ijk[I] == -1)
        || (state.ijk[J] == -1)
        || (state.ijk[K] == -1)
        || (state.ijk[I] == d_N[I])
        || (state.ijk[J] == d_N[J])
        || (state.ijk[K] == d_N[K]))
    {
        return 1;
    }
    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reflect a particle.
 */
bool Geometry::reflect(Geometry_State& state) const
{
    using def::X; using def::Y; using def::Z;

    double n0 = 0.0;
    double n1 = 0.0;
    double n2 = 0.0;
    if (state.next_ijk[X] == -1)
    {
        n0 = -1.0;
    }
    else if (state.next_ijk[X] == d_N[X])
    {
        n0 = 1.0;
    }
    else if (state.next_ijk[Y] == -1)
    {
        n1 = -1.0;
    }
    else if (state.next_ijk[Y] == d_N[Y])
    {
        n1 = 1.0;
    }
    else if (state.next_ijk[Z] == -1)
    {
        n2 = -1.0;
    }
    else if (state.next_ijk[Z] == d_N[Z])
    {
        n2 = 1.0;
    }


    // calculate the dot-product of the incoming angle and outward normal
    double dot = state.dir[X]*n0 + state.dir[Y]*n1 +
                 state.dir[Z]*n2;

    // if the dot-product != 0 then calculate the reflected angle
    if (dot != 0.0)
    {
        state.dir[X] -= 2.0 * n0 * dot;
        state.dir[Y] -= 2.0 * n1 * dot;
        state.dir[Z] -= 2.0 * n2 * dot;

        return true;
    }
    return false;
}

} // end namespace acc

//---------------------------------------------------------------------------//
//                 end of Geometry.cc
//---------------------------------------------------------------------------//
