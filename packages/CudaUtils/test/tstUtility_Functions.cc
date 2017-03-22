//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/tstUtility_Functions.cc
 * \author Tom Evans
 * \date   Mon Dec 05 14:42:05 2016
 * \brief  Tests for class Utility_Functions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Utility_Functions_Tester.hh"

//---------------------------------------------------------------------------//
// THREAD ID TESTS
//---------------------------------------------------------------------------//

TEST_F(ThreadID_Test, 1Dgrid_1Dblocks)
{
    test_1D_1D();
}

TEST_F(ThreadID_Test, 1Dgrid_2Dblocks)
{
    test_1D_2D();
}

TEST_F(ThreadID_Test, 1Dgrid_3Dblocks)
{
    test_1D_3D();
}

TEST_F(ThreadID_Test, 2Dgrid_1Dblocks)
{
    test_2D_1D();
}

TEST_F(ThreadID_Test, 2Dgrid_2Dblocks)
{
    test_2D_2D();
}

TEST_F(ThreadID_Test, 2Dgrid_3Dblocks)
{
    test_2D_3D();
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/tstUtility_Functions.cc
//---------------------------------------------------------------------------//
