//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   acc/test/PhysicsTest.hh
 * \author Seth R Johnson
 * \date   Thu Oct 30 09:35:40 2014
 * \brief  PhysicsTest class declaration.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_test_PhysicsTest_hh
#define acc_test_PhysicsTest_hh

#include "../Physics.hh"

void calc_total(
        const acc::Physics&  physics,
        const int*           matids,
        const int*           groups,
        int                  size,
        double*              total);

#endif // acc_test_PhysicsTest_hh

//---------------------------------------------------------------------------//
// end of acc/test/PhysicsTest.hh
//---------------------------------------------------------------------------//
