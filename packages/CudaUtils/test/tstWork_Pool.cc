//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/tstWork_Pool.cc
 * \author Steven Hamilton
 * \date   Wed Mar 08 13:59:17 2017
 * \brief  Tests for class Work_Pool.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include "Work_Pool_Tester.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(WorkPoolTest, test_pool)
{
    Work_Pool_Tester::test_pool();
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/tstWork_Pool.cc
//---------------------------------------------------------------------------//
