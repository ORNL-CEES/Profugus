//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/tstRTK_Cell.cc
 * \author Tom Evans
 * \date   Tue Nov 29 17:09:01 2016
 * \brief  Tests for class RTK_Cell.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RTK_Cell.hh"

#include "Nemesis/gtest/nemesis_gtest.hh"

using cuda_profugus::RTK_Cell;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RTKCellTest : public ::nemesis::Test
{
  protected:
    // >>> TYPEDEFS

  protected:
    void SetUp()
    {
        // inp_filename = nemesis::test_data_path("cuda_rtk", "inp.h5");
        // out_filename = make_unique_filename(".h5")
    }

  protected:
    // >>> DATA
    // std::string inp_filename;
    // std::string out_filename;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RTKCellTest, subtest_description)
{
    // if (nodes != 4)
    //     SKIP_TEST("This test requires 4 nodes.");
    //
    // EXPECT_EQ( 1, q.size());
    // EXPECT_SOFTEQ(1.23, q[0], 1.e-12);
    // EXPECT_SOFT_EQ(1.23, q[0]);
    // EXPECT_LT(20, q.size() + 22) << "descriptive failure message";
    // ASSERT_NE( 3, q.size()); // this will abort only the current
                                // TEST_F function

    // Test all elements in two contiguous containers
    // EXPECT_VEC_EQ(expected, actual);
    // EXPECT_VEC_SOFTEQ(expected_floats, actual_floats);

    // Print elements of a vector
    // PRINT_EXPECTED(expected)

    // for more details, see:
    // http://code.google.com/p/googletest/wiki/Primer
    // http://code.google.com/p/googletest/wiki/AdvancedGuide
    // and the Exnihilo unit test documentation
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/tstRTK_Cell.cc
//---------------------------------------------------------------------------//
