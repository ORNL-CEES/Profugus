//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstKeff_Tally_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Keff_Tally
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Keff_Tally_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(KeffTally, onegroup)
{
    Keff_Tally_Tester::test_tally(1);
}

TEST(KeffTally, twogroup)
{
    Keff_Tally_Tester::test_tally(2);
}

//---------------------------------------------------------------------------//
//                 end of tstKeff_Tally_cuda.cc
//---------------------------------------------------------------------------//
