//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/tstRTK_Geometry_cuda.cc
 * \author Tom Evans
 * \date   Fri Feb 03 09:51:07 2017
 * \brief  Tests for class RTK_Geometry_cuda.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Utils/gtest/utils_gtest.hh"

#include "RTK_Geometry_Tester.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Core, heuristic)
{
    heuristic();
}

//---------------------------------------------------------------------------//

TEST_F(Core, reflecting)
{
    reflecting();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/tstRTK_Geometry_cuda.cc
//---------------------------------------------------------------------------//
