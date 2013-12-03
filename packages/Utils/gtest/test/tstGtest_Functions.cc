//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   gtest/test/tstGtest_Functions.cc
 * \author Seth R Johnson
 * \date   Thu Apr 11 10:03:49 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/nemesis_gtest.hh"

#include "../Gtest_Functions.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(skipmessage, does_nothing)
{
    // shouldn't throw or do anything weird
    nemesis::print_skip_message();
}

//---------------------------------------------------------------------------//
TEST(IsSoftEquiv, successes)
{
    using nemesis::IsSoftEquiv;

    EXPECT_TRUE(IsSoftEquiv(
                "1.", "1. + 1.e-13", "1.e-12",
                 1. ,  1. + 1.e-13 ,  1.e-12 ));
    EXPECT_TRUE(IsSoftEquiv(
                "0.", "1.e-13", "1.e-12",
                 0. ,  1.e-13 ,  1.e-12 ));
    EXPECT_TRUE(IsSoftEquiv(
                "10.", "10.05", "1.e-2",
                 10. ,  10.05 ,  1.e-2 ));

    // no default tolerance
    EXPECT_TRUE(IsSoftEquiv(
                "1.", "1. + 1.e-13",
                 1. ,  1. + 1.e-13  ));
}

TEST(IsSoftEquiv, failures)
{
    using nemesis::IsSoftEquiv;

    EXPECT_FALSE(IsSoftEquiv(
                "1.", "1. + 1.e-13", "1.e-14",
                 1. ,  1. + 1.e-13 ,  1.e-14 ));
    EXPECT_FALSE(IsSoftEquiv(
                "0.", "2.e-13", "1.e-13",
                 0. ,  2.e-13 ,  1.e-13 ));
    EXPECT_FALSE(IsSoftEquiv(
                "10.", "10.05", "1.e-3",
                 10. ,  10.05 ,  1.e-3 ));

    // no default tolerance
    EXPECT_FALSE(IsSoftEquiv(
                "1.", "1. + 1.e-11",
                 1. ,  1. + 1.e-11  ));
}

//---------------------------------------------------------------------------//
TEST(IsVecSoftEquiv, successes)
{
    using nemesis::IsVecSoftEquiv;

    typedef std::vector<double> Vec_Dbl;

    Vec_Dbl expected, actual;
    double expected_array[] = {1.005, 0.009};

    // empty vectors
    EXPECT_TRUE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // keep the same size
    expected.push_back(1.); actual.push_back(1.005);
    EXPECT_TRUE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // test comparison against zero
    expected.push_back(0.); actual.push_back(0.009);
    EXPECT_TRUE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    EXPECT_TRUE(IsVecSoftEquiv(
                "expected_array", "actual", "0.01",
                 expected_array ,  actual ,  0.01 ));

    // no default tolerance
    expected.clear(); expected.push_back(8.);
    actual.clear(); actual.push_back(8. * (1 + 1.e-13));
    EXPECT_VEC_SOFT_EQ(expected, actual);
}

// Note: to test what the output looks like, just change EXPECT_FALSE to
// EXPECT_TRUE, or do cout << IsVecSoftEquiv(...).message() << endl;
TEST(IsVecSoftEquiv, failures)
{
    using nemesis::IsVecSoftEquiv;

    typedef std::vector<double> Vec_Dbl;

    Vec_Dbl expected, actual;

    // wrong sized vectors
    expected.push_back(1.);
    EXPECT_FALSE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // single wrong value
    actual.push_back(0.);
    EXPECT_FALSE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // truncating long expression values
    EXPECT_FALSE(IsVecSoftEquiv(
                "expectedddddddddd", "actualllllllllll", "0.01",
                 expected ,  actual ,  0.01 ));

    // multiple wrong values
    expected.push_back(100.); actual.push_back(102.);
    EXPECT_FALSE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // Ten wrong values
    expected.assign(10, 5.);
    actual.assign(10, 6.);
    EXPECT_FALSE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // Uncomment this line to test output
    //EXPECT_VEC_SOFTEQ(expected, actual, 0.01);

    // 100 wrong values (should truncate)
    expected.assign(100, 5.);
    actual.assign(100, 6.);
    EXPECT_FALSE(IsVecSoftEquiv(
                "expected", "actual", "0.01",
                 expected ,  actual ,  0.01 ));

    // A couple of wrong values in a large array
    expected.assign(200, 10.);
    actual.assign(200, 10.);
    actual[0] = 11.;
    actual[5] = 1. + 1.e-11;
    actual[100] = 3.;
    actual[150] = 2.;
    EXPECT_FALSE(IsVecSoftEquiv(
                "expected", "actual", "1.e-12",
                 expected ,  actual ,  1.e-12 ));
}

//---------------------------------------------------------------------------//
TEST(IsVecEq, successes)
{
    using nemesis::IsVecEq;

    typedef std::vector<int> Vec_Int;

    Vec_Int expected, actual;

    // empty vectors
    EXPECT_TRUE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));

    // keep the same size
    expected.push_back(1); actual.push_back(1);
    EXPECT_TRUE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));

    // test comparison against zero
    expected.push_back(0); actual.push_back(0);
    EXPECT_TRUE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));

    expected.clear(); expected.push_back(8);
    actual.clear(); actual.push_back(8);
    EXPECT_VEC_EQ(expected, actual);

    int static_array[] = {8, 0, 3};
    actual.clear();
    actual.push_back(8);
    actual.push_back(0);
    actual.push_back(3);

    EXPECT_VEC_EQ(static_array, actual);
}

// Note: to test what the output looks like, just change EXPECT_FALSE to
// EXPECT_TRUE, or do cout << IsVecEq(...).message() << endl;
TEST(IsVecEq, failures)
{
    using nemesis::IsVecEq;

    typedef std::vector<int> Vec_Int;

    Vec_Int expected, actual;

    int static_array[] = {1, 100};

    // wrong sized vectors
    expected.push_back(1);
    EXPECT_FALSE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));
    EXPECT_FALSE(IsVecEq(
                "static_array", "actual",
                 static_array ,  actual ));

    // single wrong value
    actual.push_back(0);
    EXPECT_FALSE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));

    // truncating long expression values
    EXPECT_FALSE(IsVecEq(
                "expectedddddddddd", "actualllllllllll",
                 expected ,  actual  ));

    // multiple wrong values
    expected.push_back(100.); actual.push_back(102.);
    EXPECT_FALSE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));
    EXPECT_FALSE(IsVecEq(
                "static_array", "actual",
                 static_array ,  actual ));

    // Ten wrong values
    expected.assign(10, 5.);
    actual.assign(10, 6.);
    EXPECT_FALSE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));

    // Uncomment this line to test output
    //EXPECT_VEC_EQ(expected, actual);

    // 100 wrong values (should truncate)
    expected.assign(100, 5.);
    actual.assign(100, 6.);
    EXPECT_FALSE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));

    // A couple of wrong values in a large array
    expected.assign(200, 10);
    actual.assign(200, 10);
    actual[0] = 1;
    actual[5] = 1;
    actual[100] = 3;
    actual[150] = 2;
    EXPECT_FALSE(IsVecEq(
                "expected", "actual",
                 expected ,  actual ));
}

//---------------------------------------------------------------------------//
//                        end of tstGtest_Functions.cc
//---------------------------------------------------------------------------//
