//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/test/PhysicsTest.cc
 * \author Seth R Johnson
 * \date   Thu Oct 30 09:35:40 2014
 * \brief  PhysicsTest class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PhysicsTest.hh"

void calc_total(
        const acc::Physics& physics,
        const int*     matids,
        const int*     groups,
        int            size,
        double*        total)
{
#pragma acc parallel loop \
        present(physics) \
        copyin(groups[0:size], matids[0:size]) \
        copyout(total[0:size])
    for (int i = 0; i < size; ++i)
    {
        total[i] = physics.total(matids[i], groups[i]);
    }
}

//---------------------------------------------------------------------------//
// end of acc/test/PhysicsTest.cc
//---------------------------------------------------------------------------//
