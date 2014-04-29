//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Functions.cc
 * \author Seth R Johnson
 * \date   Fri Jan 18 19:18:14 2013
 * \brief  RTK_Functions implementation
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Functions.hh"

#include "harness/Soft_Equivalence.hh"
#include "Definitions.hh"

using def::Space_Vector;

namespace profugus
{

/*!
 * \brief Move a particle being initialized to its intersection with the Array
 *
 * This function uses the state's current d_r and d_dir, transporting it to the
 * boundaries of this array object. If the particle never intersects, the
 * position is set to an outer edge, and the exiting face flags are set.
 *
 * \param[in] lower  Lower coordinate of the array/pin/etc.
 * \param[in] upper  Upper coordinate of the array/pin/etc.
 * \param[out] state Geometry state
 */
void move_from_outside(
        const def::Space_Vector& lower,
        const def::Space_Vector& upper,
        RTK_State&               state)
{
    const Space_Vector& omega = state.d_dir;
    Space_Vector& r           = state.d_r;
    using def::X; using def::Y; using def::Z;

    Require(nemesis::soft_equiv(vector_magnitude(omega), 1., 1.e-6));
    Require(   (omega[X] < 0 ? r[X] > lower[X] : r[X] < upper[X])
            || (omega[Y] < 0 ? r[Y] > lower[Y] : r[Y] < upper[Y])
            || (omega[Z] < 0 ? r[Z] > lower[Z] : r[Z] < upper[Z]));

    // We want to transport to the furthest possible face, checking only faces
    // that are intersectable by the particle and closest to the particle
    double temp_distance;
    double distance = 0.;

    // Check X faces if outside the extents
    if      (omega[X] < 0.0 && r[X] > upper[X])
    {
        temp_distance = (upper[X] - r[X]) / omega[X];
        if (temp_distance > distance)
            distance = temp_distance;
    }
    else if (omega[X] > 0.0 && r[X] < lower[X])
    {
        temp_distance = (lower[X] - r[X]) / omega[X];
        if (temp_distance > distance)
            distance = temp_distance;
    }
    // Check Y faces if outside the extents
    if      (omega[Y] < 0.0 && r[Y] > upper[Y])
    {
        temp_distance = (upper[Y] - r[Y]) / omega[Y];
        if (temp_distance > distance)
            distance = temp_distance;
    }
    else if (omega[Y] > 0.0 && r[Y] < lower[Y])
    {
        temp_distance = (lower[Y] - r[Y]) / omega[Y];
        if (temp_distance > distance)
            distance = temp_distance;
    }
    // Check Z faces if outside the extents
    if      (omega[Z] < 0.0 && r[Z] > upper[Z])
    {
        temp_distance = (upper[Z] - r[Z]) / omega[Z];
        if (temp_distance > distance)
            distance = temp_distance;
    }
    else if (omega[Z] > 0.0 && r[Z] < lower[Z])
    {
        temp_distance = (lower[Z] - r[Z]) / omega[Z];
        if (temp_distance > distance)
            distance = temp_distance;
    }

    // Distance might be zero if the particle is outside the problem and
    // parallel to the face
    Check(distance > 0. || omega[X] == 0. || omega[Y] == 0. || omega[Z] == 0.);

    // Now move the particle, with a little extra push to make sure it's just
    // inside the boundary
    r[X] += distance * omega[X];
    r[Y] += distance * omega[Y];
    r[Z] += distance * omega[Z];

    // Now again, check if it missed our outer boundaries and move it to inside
    // the array
    // Check(state.escaping_face == RTK_State::NONE);
    if      (omega[X] >= 0.0 && r[X] > upper[X])
    {
        state.escaping_face = RTK_State::PLUS_X;
        r[X] = upper[X];
    }
    else if (omega[X] <= 0.0 && r[X] < lower[X])
    {
        state.escaping_face = RTK_State::MINUS_X;
        r[X] = lower[X];
    }

    if      (omega[Y] >= 0.0 && r[Y] > upper[Y])
    {
        state.escaping_face = RTK_State::PLUS_Y;
        r[Y] = upper[Y];
    }
    else if (omega[Y] <= 0.0 && r[Y] < lower[Y])
    {
        state.escaping_face = RTK_State::MINUS_Y;
        r[Y] = lower[Y];
    }

    if      (omega[Z] >= 0.0 && r[Z] > upper[Z])
    {
        state.escaping_face = RTK_State::PLUS_Z;
        r[Z] = upper[Z];
    }
    else if (omega[Z] <= 0.0 && r[Z] < lower[Z])
    {
        state.escaping_face = RTK_State::MINUS_Z;
        r[Z] = lower[Z];
    }

    Ensure((state.d_r[X] >= lower[X]) && (state.d_r[X] <= upper[X]));
    Ensure((state.d_r[Y] >= lower[Y]) && (state.d_r[Y] <= upper[Y]));
    Ensure((state.d_r[Z] >= lower[Z]) && (state.d_r[Z] <= upper[Z]));
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RTK_Functions.cc
//---------------------------------------------------------------------------//
