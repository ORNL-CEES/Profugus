//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/Change_Direction.acc.cc
 * \author Seth R Johnson
 * \date   Thu Oct 30 13:38:47 2014
 * \brief  Change_Direction.acc class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Change_Direction.hh"

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <cmath>
#endif

#include "Geometry_State.hh"
#include "core/mc/Definitions.hh"
#include "core/geometry/Definitions.hh"

namespace acc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Change the direction through an angle.
 */
void change_direction(
        double          costheta,
        double          phi,
        Geometry_State& state)
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
} // end namespace acc

//---------------------------------------------------------------------------//
// end of acc/Change_Direction.acc.cc
//---------------------------------------------------------------------------//
