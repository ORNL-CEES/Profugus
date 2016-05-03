//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstBox_Shape_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 15:22:15 2016
 * \brief  Test for Box_Shape
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Utils/gtest/utils_gtest.hh"
#include "Box_Shape_Tester.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(BoxShape, inside)
{
    Box_Shape_Tester::test_inside();
}

TEST(BoxShape, sample)
{
    Box_Shape_Tester::test_sample();
}

//---------------------------------------------------------------------------//
//                 end of tstBox_Shape.cc
//---------------------------------------------------------------------------//
