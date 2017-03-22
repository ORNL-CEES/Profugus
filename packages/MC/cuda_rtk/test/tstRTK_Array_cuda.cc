//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/tstRTK_Array_cuda.cc
 * \author Tom Evans
 * \date   Wed Jan 04 23:29:04 2017
 * \brief  Tests for class RTK_Array_cuda.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Utils/gtest/utils_gtest.hh"

#include "RTK_Array_Tester.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// These tests are all derived from tstRTK_Array.cc

TEST_F(SimpleLattice, init)
{
}

//---------------------------------------------------------------------------//

TEST_F(SimpleCore, init)
{
}

//---------------------------------------------------------------------------//

TEST_F(RegCore, track)
{
    core->complete(0.0, 0.0, 0.0);
    run_test();
}

//---------------------------------------------------------------------------//

TEST_F(RegCore, reflect)
{
    Core::Vec_Int reflect = {1, 0, 1, 0, 1, 0};
    core->set_reflecting(reflect);
    core->complete(0.0, 0.0, 0.0);
    reflect_test();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/tstRTK_Array_cuda.cc
//---------------------------------------------------------------------------//
