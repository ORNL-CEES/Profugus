//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   acc/Change_Direction.hh
 * \author Seth R Johnson
 * \date   Thu Oct 30 13:25:59 2014
 * \brief  Change_Direction class declaration.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_Change_Direction_hh
#define acc_Change_Direction_hh

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <cmath>
#endif

#include "Geometry_State.hh"

namespace acc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Change the direction through an angle.
 */
#pragma acc routine seq
void change_direction(
        double          costheta,
        double          phi,
        Geometry_State& state);

//---------------------------------------------------------------------------//
} // end namespace acc

#endif // acc_Change_Direction_hh

//---------------------------------------------------------------------------//
// end of acc/Change_Direction.hh
//---------------------------------------------------------------------------//
