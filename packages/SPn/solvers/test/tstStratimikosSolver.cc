//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstStratimikosSolver.cc
 * \author 9te
 * \date   Mon Feb 17 21:12:05 2014
 * \brief
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/nemesis_gtest.hh"

#include "../StratimikosSolver.hh"

using ::StratimikosSolver;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class StratimikosSolverTest : public ::testing::Test
{
  protected:
    // Typedefs usable inside the test fixture

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> Data that get re-initialized between tests
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(StratimikosSolverTest, subtest_description)
{
    // EXPECT_EQ( 1, q.size());
    // EXPECT_SOFTEQ(1.23, q[0], 1.e-12);
    // EXPECT_SOFT_EQ(1.23, q[0]);
    // EXPECT_LT(20, q.size() + 22) << "descriptive failure message";
    // ASSERT_NE( 3, q.size()); // this will abort only the current
                                // TEST_F function

    // Test all elements in two contiguous containers
    // EXPECT_VEC_EQ(expected, actual);
    // EXPECT_VEC_SOFTEQ(expected_floats, actual_floats);

    // for more details, see:
    // http://code.google.com/p/googletest/wiki/Primer
    // http://code.google.com/p/googletest/wiki/AdvancedGuide
}

//---------------------------------------------------------------------------//

TEST_F(StratimikosSolverTest, another_subtest)
{
    /* * */
}

//---------------------------------------------------------------------------//
//                 end of tstStratimikosSolver.cc
//---------------------------------------------------------------------------//
