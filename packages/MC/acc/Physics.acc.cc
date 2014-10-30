//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/Physics.acc.cc
 * \author Seth R Johnson
 * \date   Wed Oct 29 10:32:32 2014
 * \brief  Physics class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics.hh"

#include "harness/DBC.hh"

namespace acc
{
//---------------------------------------------------------------------------//
void Physics::complete()
{
    REQUIRE(d_num_mats > 0);
    REQUIRE(d_num_groups > 0);
    REQUIRE(d_total != 0);
    REQUIRE(d_nusigf != 0);
    REQUIRE(d_scatter != 0);
    REQUIRE(d_scatter_ratio != 0);
    REQUIRE(d_fissionable != 0);
    int ne = dv_total.size();
    int ns = dv_scatter.size();
#pragma acc enter data \
    copyin(this)
#pragma acc enter data \
    copyin(d_total[0:ne], d_nusigf[0:ne], d_scatter[0:ns],\
           d_scatter_ratio[0:ne], d_fissionable[0:d_num_mats])
}

//---------------------------------------------------------------------------//
Physics::~Physics()
{
#pragma acc exit data \
    delete(d_total, d_nusigf, d_scatter, d_scatter_ratio, d_fissionable)
#pragma acc exit data \
    delete(this)
}

//---------------------------------------------------------------------------//
} // end namespace acc

//---------------------------------------------------------------------------//
// end of acc/Physics.acc.cc
//---------------------------------------------------------------------------//
