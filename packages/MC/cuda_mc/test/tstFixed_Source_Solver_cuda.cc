//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstFixed_Solver_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Fixed_Solver
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fixed_Solver_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(FixedSolver, five_group)
{
    Fixed_Solver_Tester::test_transport(5);
}

TEST(FixedSolver, three_group)
{
    Fixed_Solver_Tester::test_transport(3);
}

TEST(FixedSolver, one_group)
{
    Fixed_Solver_Tester::test_transport(1);
}

//---------------------------------------------------------------------------//
//                 end of tstFixed_Solver_cuda.cc
//---------------------------------------------------------------------------//
