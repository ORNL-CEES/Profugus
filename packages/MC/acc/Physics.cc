//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/Physics.cc
 * \author Seth R Johnson
 * \date   Wed Oct 29 10:32:32 2014
 * \brief  Physics class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics.hh"

namespace acc
{
//---------------------------------------------------------------------------//
void Physics::complete()
{
    Require(d_num_mats > 0);
    Require(d_num_groups > 0);
    Require(d_total != 0);
    Require(d_nusigf != 0);
    Require(d_scatter != 0);
    int ne = dv_total.size();
    int ns = dv_scatter.size();
#pragma acc enter data copyin(this)
#pragma acc enter data copyin(d_total[0:ne], d_nusigf[0:ne], d_scatter[0:ns)
}
//---------------------------------------------------------------------------//
} // end namespace acc

//---------------------------------------------------------------------------//
// end of acc/Physics.cc
//---------------------------------------------------------------------------//
