//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/tstRTK_State_cuda.cc
 * \author Tom Evans
 * \date   Fri Nov 18 15:25:20 2016
 * \brief  Tests for class RTK_State.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TEST PROTOTYPES
// see RTK_State_Tester.cu for implementations
//---------------------------------------------------------------------------//

namespace rtk_state_test
{

void test_on_device();
void test_to_device();

} // end namespace rtk_state_test

//---------------------------------------------------------------------------//

TEST(RTK_StateTest, on_device)
{
    rtk_state_test::test_on_device();
}

//---------------------------------------------------------------------------//

TEST(RTK_StateTest, to_device)
{
    rtk_state_test::test_to_device();
}

//---------------------------------------------------------------------------//
// end of tstRTK_State_cuda.cc
//---------------------------------------------------------------------------//
