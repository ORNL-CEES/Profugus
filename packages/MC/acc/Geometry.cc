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

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <cmath>
#endif

#include "Geometry.hh"

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

#ifdef _OPENACC
    state.ijk[0] = r[0] / (d_x[1] - d_x[0]);
    state.ijk[1] = r[1] / (d_y[1] - d_y[0]);
    state.ijk[2] = r[2] / (d_z[1] - d_z[0]);
#else
    state.ijk[0] = lower_bound(d_x, d_x + d_N[0] + 1, state.pos[0]) - d_x - 1;
    state.ijk[1] = lower_bound(d_y, d_y + d_N[1] + 1, state.pos[1]) - d_y - 1;
    state.ijk[2] = lower_bound(d_z, d_z + d_N[2] + 1, state.pos[2]) - d_z - 1;
#endif
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
void Geometry::change_direction(double          costheta,
                                double          phi,
                                Geometry_State& state) const
{
#ifndef _OPENACC
    using std::cos; using std::sin; using std::sqrt;
#endif
    using def::I; using def::J; using def::K;

     // cos/sin factors
    const double cosphi   = cos(phi);
    const double sinphi   = sin(phi);
    const double sintheta = sqrt(1.0 - costheta * costheta);

    // make a copy of the d_work direction
    double o0 = state.dir[0];
    double o1 = state.dir[1];
    double o2 = state.dir[2];

    // calculate alpha
    const double alpha = sqrt(1.0 - o2 * o2);

    // now transform into new cooordinate direction; degenerate case first
    if (alpha < 1.e-6)
    {
        state.dir[I] = sintheta * cosphi;
        state.dir[J] = sintheta * sinphi;
        state.dir[K] = (o2 < 0.0 ? -1.0 : 1.0) * costheta;
    }

    // do standard transformation
    else
    {
        // calculate inverse of alpha
        const double inv_alpha = 1.0 / alpha;

        // calculate new z-direction
        state.dir[K] = o2 * costheta - alpha * sintheta * cosphi;

        // calculate new x-direction
        state.dir[I] = o0 * costheta + inv_alpha * (
            o0 * o2 * sintheta * cosphi - o1 * sintheta * sinphi);

        // calculate new y-direction
        state.dir[J] = o1 * costheta + inv_alpha * (
            o1 * o2 * sintheta * cosphi + o0 * sintheta * sinphi);
    }

    // normalize the particle to avoid roundoff errors
    double norm = 1.0 / sqrt(state.dir[I] * state.dir[I] +
                             state.dir[J] * state.dir[J] +
                             state.dir[K] * state.dir[K]);
    state.dir[0] *= norm;
    state.dir[1] *= norm;
    state.dir[2] *= norm;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the boundary state.
 */
int Geometry::boundary_state(const Geometry_State &state) const
{
    using def::I; using def::J; using def::K;

    for (int d = 0; d < 3; ++d)
    {
        if (d_b[d*2] && state.next_ijk[d] == -1)
        {
            return profugus::geometry::REFLECT;
        }
        else if (d_b[d*2+1] && state.next_ijk[d] == d_N[d])
        {
            return profugus::geometry::REFLECT;
        }
    }

    if ((state.ijk[I] == -1)
        || (state.ijk[J] == -1)
        || (state.ijk[K] == -1)
        || (state.ijk[I] == d_N[I])
        || (state.ijk[J] == d_N[J])
        || (state.ijk[K] == d_N[K]))
    {
        return profugus::geometry::OUTSIDE;
    }
    return profugus::geometry::INSIDE;
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
