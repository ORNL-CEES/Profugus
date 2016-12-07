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

TEST_F(Single_Shell, construction)
{
    construct();
}

//---------------------------------------------------------------------------//

TEST_F(Single_Shell, tracking)
{
    track();
}

//---------------------------------------------------------------------------//

TEST_F(Multi_Shell, construction)
{
    construct();
}

//---------------------------------------------------------------------------//

TEST_F(Multi_Shell, tracking)
{
    track();
}

//---------------------------------------------------------------------------//

TEST_F(Multi_Shell, multiseg_construction)
{
    multiseg_construct();
}

//---------------------------------------------------------------------------//

TEST_F(Multi_Shell, multisegment_tracking)
{
    multiseg_track();
}

//---------------------------------------------------------------------------//

TEST_F(Empty, square)
{
    square();
}

//---------------------------------------------------------------------------//

TEST_F(Empty, rectangle)
{
    rectangle();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/tstRTK_Cell.cc
//---------------------------------------------------------------------------//
