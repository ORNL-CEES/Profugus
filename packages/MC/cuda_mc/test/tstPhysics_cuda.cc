//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstPhysics_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Physics
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Physics, total)
{
    Physics_Tester::test_total();
}

TEST(Physics, collide)
{
    Physics_Tester::test_collide();
}

//---------------------------------------------------------------------------//
//                 end of tstPhysics_cuda.cc
//---------------------------------------------------------------------------//
