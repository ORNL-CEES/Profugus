//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstSoft_Comparators.cc
 * \author Gregory G. Davidson
 * \date   Wed Sep 11 21:33:15 2013
 * \brief  Tests the Soft_Comparators
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Soft_Comparators.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class Soft_ComparatorsTest : public ::testing::Test
{
  protected:
    // Typedefs usable inside the test fixture

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        epsilon = 0.01;
        val_1 = 1.0;
        val_2 = 1.0001;
        val_3 = 1.02;
    }

  protected:
    // >>> Data that get re-initialized between tests
    double epsilon;
    double val_1;
    double val_2;
    double val_3;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Soft_ComparatorsTest, subtest_description)
{
    EXPECT_TRUE(  profugus::soft_is_equal(val_1, val_2, epsilon) );
    EXPECT_TRUE( !profugus::soft_is_equal(val_1, val_3, epsilon) );

    // Test IsNotEqual
    EXPECT_TRUE(  profugus::soft_is_not_equal(val_1, val_3, epsilon) );
    EXPECT_TRUE( !profugus::soft_is_not_equal(val_1, val_2, epsilon) );

    // Test IsLess
    EXPECT_TRUE(  profugus::soft_is_less(val_1, val_3, epsilon) );
    EXPECT_TRUE( !profugus::soft_is_less(val_1, val_2, epsilon) );

    // Test IsLessEqual
    EXPECT_TRUE( profugus::soft_is_less_equal(val_1, val_2, epsilon) );
    EXPECT_TRUE( profugus::soft_is_less_equal(val_1, val_3, epsilon) );

    // Test IsGreater
    EXPECT_TRUE(  profugus::soft_is_greater(val_3, val_1, epsilon) );
    EXPECT_TRUE( !profugus::soft_is_greater(val_2, val_1, epsilon) );

    // Test IsGreaterEqual
    EXPECT_TRUE( profugus::soft_is_greater_equal(val_2, val_1, epsilon) );
    EXPECT_TRUE( profugus::soft_is_greater_equal(val_3, val_1, epsilon) );

    // Test is_within
    EXPECT_TRUE(  profugus::soft_is_within(std::make_pair(1.0, 2.0),
                                           std::make_pair(0.99, 2.01), 0.01) );
    EXPECT_TRUE( !profugus::soft_is_within(std::make_pair(1.0, 2.0),
                                           std::make_pair(1.2, 2.0), 0.01) );
    EXPECT_TRUE( !profugus::soft_is_within(std::make_pair(1.0, 2.0),
                                           std::make_pair(1.0, 1.8), 0.01) );
}

//---------------------------------------------------------------------------//
//                        end of tstSoft_Comparators.cc
//---------------------------------------------------------------------------//
