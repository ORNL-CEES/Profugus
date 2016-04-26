//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   gtest/Gtest_Functions.hh
 * \author Seth R Johnson
 * \date   Tue Apr 02 11:41:41 2013
 * \brief  Gtest_Functions class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef gtest_Gtest_Functions_hh
#define gtest_Gtest_Functions_hh

#include "gtest.h"

namespace profugus
{
//===========================================================================//

// Implementation for the harness' "main" function.
int gtest_main(int argc, char *argv[]);

// Print the "skip" message from the skip macro
void print_skip_message();

//---------------------------------------------------------------------------//
// Custom error mesages for relative error soft equiavelence
::testing::AssertionResult IsSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* eps_expr,
        double expected,
        double actual,
        double eps);

//! Soft equivalence assertion with default argument
inline ::testing::AssertionResult IsSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        double expected,
        double actual)
{
    return IsSoftEquiv(
            expected_expr, actual_expr, "1e-12",
            expected     , actual     , 1.e-12);
}

//---------------------------------------------------------------------------//
// Helper function: analyze using raw pointers for compilable flexibility
::testing::AssertionResult IsIterSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* eps_expr,
        const double* expected_begin, const double* expected_end,
        const double* actual_begin,   const double* actual_end,
        double eps);

//---------------------------------------------------------------------------//
// Custom vector soft equivalence comparison
template<class Container_E, class Container_A>
::testing::AssertionResult IsVecSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* eps_expr,
        const Container_E& expected,
        const Container_A& actual,
        double eps)
{
    // Get pointers using bracket operator
    const double* expected_begin = &expected[0];
    const double* expected_end   = &expected[expected.size() - 1] + 1;
    const double* actual_begin   = &actual[0];
    const double* actual_end     = &actual[actual.size() - 1] + 1;

    // Ensure continuity of underlying data
    if (expected_end - expected_begin != static_cast<int>(expected.size()))
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Noncontiguous container: " << expected_expr;
        return failure;
    }
    if (actual_end - actual_begin != static_cast<int>(actual.size()))
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Noncontiguous container: " << actual_expr;
        return failure;
    }

    return IsIterSoftEquiv(
            expected_expr, actual_expr, eps_expr,
            expected_begin, expected_end,
            actual_begin, actual_end,
            eps);
}

//---------------------------------------------------------------------------//
// Custom vector equality comparison against a statically sized C array
template<std::size_t N, class Container_A>
inline ::testing::AssertionResult IsVecSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* eps_expr,
        const double (&expected)[N],
        const Container_A& actual,
        double eps)
{
    // Get pointers using bracket operator
    const double* expected_begin = expected;
    const double* expected_end   = expected + N;
    const double* actual_begin   = &actual[0];
    const double* actual_end     = &actual[actual.size() - 1] + 1;

    // Ensure continuity of underlying data
    if (actual_end - actual_begin != static_cast<int>(actual.size()))
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Noncontiguous container: " << actual_expr;
        return failure;
    }

    return IsIterSoftEquiv(
            expected_expr, actual_expr, eps_expr,
            expected_begin, expected_end,
            actual_begin, actual_end,
            eps);
}

//---------------------------------------------------------------------------//
// Helper function: analyze using raw pointers for compilable flexibility
template<class T>
::testing::AssertionResult IsIterEq(
        const char* expected_expr,
        const char* actual_expr,
        const T* expected, const T* expected_end,
        const T* actual,   const T* actual_end);

//---------------------------------------------------------------------------//
// Custom vector equality comparison
template<class Container_E, class Container_A>
inline ::testing::AssertionResult IsVecEq(
        const char* expected_expr,
        const char* actual_expr,
        const Container_E& expected,
        const Container_A& actual)
{
    typedef const typename Container_E::value_type* const_ptr_E;
    typedef const typename Container_A::value_type* const_ptr_A;

    // Get pointers using bracket operator
    const_ptr_E expected_begin = &expected[0];
    const_ptr_E expected_end   = &expected[expected.size() - 1] + 1;
    const_ptr_A actual_begin   = &actual[0];
    const_ptr_A actual_end     = &actual[actual.size() - 1] + 1;

    // Ensure continuity of underlying data
    if (expected_end - expected_begin != static_cast<int>(expected.size()))
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Noncontiguous container: " << expected_expr;
        return failure;
    }
    if (actual_end - actual_begin != static_cast<int>(actual.size()))
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Noncontiguous container: " << actual_expr;
        return failure;
    }

    return IsIterEq(
            expected_expr, actual_expr,
            expected_begin, expected_end,
            actual_begin, actual_end);
}

//---------------------------------------------------------------------------//
// Custom vector equality comparison against a statically sized C array
template<typename T, std::size_t N, class Container_A>
inline ::testing::AssertionResult IsVecEq(
        const char* expected_expr,
        const char* actual_expr,
        const T (&expected)[N],
        const Container_A& actual)
{
    typedef const T* const_ptr_E;
    typedef const typename Container_A::value_type* const_ptr_A;

    // Get pointers using bracket operator
    const_ptr_E expected_begin = expected;
    const_ptr_E expected_end   = expected + N;
    const_ptr_A actual_begin   = &actual[0];
    const_ptr_A actual_end     = &actual[actual.size() - 1] + 1;

    // Ensure continuity of underlying data
    if (actual_end - actual_begin != static_cast<int>(actual.size()))
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Noncontiguous container: " << actual_expr;
        return failure;
    }

    return IsIterEq(
            expected_expr, actual_expr,
            expected_begin, expected_end,
            actual_begin, actual_end);
}

//===========================================================================//
} // end namespace profugus

//---------------------------------------------------------------------------//
// TESTING MACROS
//---------------------------------------------------------------------------//

//! Soft equivalence macro
#define EXPECT_SOFTEQ(expected, actual, rel_error) \
    EXPECT_PRED_FORMAT3(::profugus::IsSoftEquiv, expected, actual, rel_error)

//! Soft equivalence macro with default argument (like google test)
#define EXPECT_SOFT_EQ(expected, actual) \
    EXPECT_PRED_FORMAT2(::profugus::IsSoftEquiv, expected, actual)

//! Soft equivalence macro for containers of doubles
#define EXPECT_VEC_SOFTEQ(expected, actual, rel_error) \
    EXPECT_PRED_FORMAT3(::profugus::IsVecSoftEquiv, expected, actual, rel_error)

//! Soft equivalence macro for containers of doubles
#define EXPECT_VEC_SOFT_EQ(expected, actual) \
    EXPECT_VEC_SOFTEQ(expected, actual, 1e-12)

//! Equivalence macro for containers of integers or whatever
#define EXPECT_VEC_EQ(expected, actual) \
    EXPECT_PRED_FORMAT2(::profugus::IsVecEq, expected, actual)

//! Macro to skip a unit test and print a colored warning
#define SKIP_TEST(REASON_STRING_STREAM) \
    do { \
        ::profugus::print_skip_message(); \
        std::cout << REASON_STRING_STREAM << std::endl; \
        return; \
    } while (0)


#endif // gtest_Gtest_Functions_hh

//---------------------------------------------------------------------------//
//              end of gtest/Gtest_Functions.hh
//---------------------------------------------------------------------------//
