//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/geometry/RTK_Functions.hh
 * \author Seth R Johnson
 * \date   Fri Jan 18 19:18:14 2013
 * \brief  RTK_Functions class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_geometry_RTK_Functions_hh
#define MC_geometry_RTK_Functions_hh

#include "utils/Definitions.hh"
#include "RTK_State.hh"

namespace profugus
{

// Return the sign of a value
// Taken from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// Transport to the array boundaries from a point outside the geometry.
void move_from_outside(
        const def::Space_Vector& lower,
        const def::Space_Vector& upper,
        RTK_State&               state);

} // end namespace profugus

#endif // MC_geometry_RTK_Functions_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_Functions.hh
//---------------------------------------------------------------------------//
