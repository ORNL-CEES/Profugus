//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstKCode_Solver_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for KCode_Solver
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KCode_Solver_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(KCodeSolver, five_group)
{
    KCode_Solver_Tester::test_transport(5);
}

TEST(KCode_Solver_cudaTest, three_group)
{
    KCode_Solver_Tester::test_transport(3);
}

TEST(KCodeSolver, one_group)
{
    KCode_Solver_Tester::test_transport(1);
}

//---------------------------------------------------------------------------//
//                 end of tstKCode_Solver_cuda.cc
//---------------------------------------------------------------------------//
