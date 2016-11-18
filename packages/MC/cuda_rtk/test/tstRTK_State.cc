//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/tstRTK_State.cc
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

void test_on_device();
void test_to_device();

//---------------------------------------------------------------------------//

TEST(RTK_StateTest, on_device)
{
    test_on_device();
}

//---------------------------------------------------------------------------//

TEST(RTK_StateTest, to_device)
{
    test_to_device();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/tstRTK_State.cc
//---------------------------------------------------------------------------//
