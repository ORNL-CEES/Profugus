//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstTypes.cc
 * \author 9te
 * \date   Tue Dec 03 21:48:11 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/nemesis_gtest.hh"

#include "../Types.hh"

using ::Types;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class TypesTest : public ::testing::Test
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

TEST_F(TypesTest, subtest_description)
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

TEST_F(TypesTest, another_subtest)
{
    /* * */
}

//---------------------------------------------------------------------------//
//                 end of tstTypes.cc
//---------------------------------------------------------------------------//
