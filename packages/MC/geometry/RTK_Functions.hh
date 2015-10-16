//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Functions.hh
 * \author Seth R Johnson
 * \date   Fri Jan 18 19:18:14 2013
 * \brief  RTK_Functions class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_RTK_Functions_hh
#define geometry_RTK_Functions_hh

#include "utils/Definitions.hh"
#include "RTK_State.hh"

namespace profugus
{

// Transport to the array boundaries from a point outside the geometry.
void move_from_outside(
        const def::Space_Vector& lower,
        const def::Space_Vector& upper,
        RTK_State&               state);

} // end namespace profugus

#endif // geometry_RTK_Functions_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_Functions.hh
//---------------------------------------------------------------------------//
