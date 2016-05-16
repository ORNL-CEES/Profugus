//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/gtest/utils_gtest.hh
 * \author Seth R Johnson
 * \date   Thu Aug 09 10:14:06 2012
 * \brief  Load google test macros
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * \warning This file should *only* be included from a unit test .cc file! It
 * provides the implementation of main().
 */
//---------------------------------------------------------------------------//

#ifndef Utils_gtest_utils_gtest_hh
#define Utils_gtest_utils_gtest_hh

#include <Utils/config.h>
#include "gtest.h"
#include "Gtest_Functions.hh"
#include "Test.hh"

#include <iostream>

#include "Utils/harness/DBC.hh"
#include "Utils/harness/Soft_Equivalence.hh"
#include "Utils/comm/global.hh"

// Import some common stuff into the global namespace
using std::cout;
using std::cerr;
using std::endl;
using profugus::soft_equiv;

//---------------------------------------------------------------------------//
// TESTING MACROS
//---------------------------------------------------------------------------//

//! Soft equivalence macro
#define EXPECT_SOFT_EQ(expected, actual) \
    EXPECT_PRED_FORMAT2(::profugus::IsSoftEquiv, expected, actual)

//! Soft equivalence macro with relative error
#define EXPECT_SOFTEQ(expected, actual, rel_error) \
    EXPECT_PRED_FORMAT3(::profugus::IsSoftEquiv, expected, actual, rel_error)

//! Strong soft equivalence macro with relative error
#define ASSERT_SOFTEQ(expected, actual, rel_error) \
    ASSERT_PRED_FORMAT3(::profugus::IsSoftEquiv, expected, actual, rel_error)

//! Soft equivalence macro with relative and absolute error
#define EXPECT_CLOSE(expected, actual, rel_error, abs_thresh) \
    EXPECT_PRED_FORMAT4(::profugus::IsSoftEquiv, \
                        expected, actual, rel_error, abs_thresh)

//! Container soft equivalence macro
#define EXPECT_VEC_SOFT_EQ(expected, actual) \
    EXPECT_PRED_FORMAT2(::profugus::IsVecSoftEquiv, expected, actual)

//! Container soft equivalence macro with relative error
#define EXPECT_VEC_SOFTEQ(expected, actual, rel_error) \
    EXPECT_PRED_FORMAT3(::profugus::IsVecSoftEquiv, expected, actual, rel_error)

//! Container soft equivalence macro with relative and absolute error
#define EXPECT_VEC_CLOSE(expected, actual, rel_error, abs_thresh) \
    EXPECT_PRED_FORMAT4(::profugus::IsVecSoftEquiv, \
                        expected, actual, rel_error, abs_thresh)

//! Container equality macro
#define EXPECT_VEC_EQ(expected, actual) \
    EXPECT_PRED_FORMAT2(::profugus::IsVecEq, expected, actual)

//! Macro to skip a unit test and print a colored warning
#define SKIP_TEST(REASON_STRING_STREAM) \
    do { \
        ::profugus::print_skip_message(); \
        std::cout << REASON_STRING_STREAM << std::endl; \
        return; \
    } while (0)

//! Print the given container as an array for regression testing
#define PRINT_EXPECTED(data) ::profugus::print_expected(data, #data)

//! Make expect true tests for smart pointers.
#define EXPECT_SP_TRUE(condition) \
    GTEST_TEST_BOOLEAN_(static_cast<bool>(condition), #condition, false, true, \
                        GTEST_NONFATAL_FAILURE_)

//! Make expect false tests for smart pointers.
#define EXPECT_SP_FALSE(condition) \
    GTEST_TEST_BOOLEAN_(!static_cast<bool>(condition), #condition, true, false,\
                        GTEST_NONFATAL_FAILURE_)

//! Make assert true tests for smart pointers.
#define ASSERT_SP_TRUE(condition) \
    GTEST_TEST_BOOLEAN_(static_cast<bool>(condition), #condition, false, true, \
                        GTEST_FATAL_FAILURE_)

//! Make assert false tests for smart pointers.
#define ASSERT_SP_FALSE(condition) \
    GTEST_TEST_BOOLEAN_(!static_cast<bool>(condition), #condition, true, false,\
                        GTEST_FATAL_FAILURE_)

//---------------------------------------------------------------------------//
// MAIN FUNCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Launch the google test harness from this unit test
 *
 * \warning We are intentionally putting a non-inlined function implementation
 * in the header! This is so that we can simply include
 * "gtest/utils_gtest.hh" and not have to rewrite this section of code.
 *
 * We have to use this method rather than linking in a separate "gtest-main"
 * class because of a linker bug in MS Visual Studio that deletes main()
 * functions in shared libraries.
 */
int main(int argc, char *argv[])
{
    profugus::gtest_main(argc, argv);
}

//---------------------------------------------------------------------------//

#endif // Utils_gtest_utils_gtest_hh

//---------------------------------------------------------------------------//
//              end of gtest/utils_gtest.hh
//---------------------------------------------------------------------------//
