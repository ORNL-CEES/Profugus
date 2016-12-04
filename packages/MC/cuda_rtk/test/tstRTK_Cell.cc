//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/tstRTK_Cell.cc
 * \author Tom Evans
 * \date   Tue Nov 29 17:09:01 2016
 * \brief  Tests for class RTK_Cell.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Utils/gtest/utils_gtest.hh"

#include "RTK_Cell_Tester.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// These tests are all derived from tstRTK_Cell.cc

TEST_F(Single_Shell, execution)
{
    run();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/tstRTK_Cell.cc
//---------------------------------------------------------------------------//
