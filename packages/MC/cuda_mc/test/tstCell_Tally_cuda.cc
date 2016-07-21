//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstCell_Tally_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Cell_Tally
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//


#include "Cell_Tally_Tester.hh"
#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(CellTally, tally_1_batch)
{
    Cell_Tally_Tester::test_tally(1);
}

TEST(CellTally, tally_10_batch)
{
    Cell_Tally_Tester::test_tally(10);
}

TEST(CellTally, tally_100_batch)
{
    Cell_Tally_Tester::test_tally(100);
}

//---------------------------------------------------------------------------//
//                 end of tstCell_Tally_cuda.cc
//---------------------------------------------------------------------------//
