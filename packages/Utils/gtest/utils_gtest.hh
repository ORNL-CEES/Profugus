//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   gtest/utils_gtest.hh
 * \author Seth R Johnson
 * \date   Thu Aug 09 10:14:06 2012
 * \brief  Load google test macros
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * \warning This file should *only* be included from a unit test .cc file! It
 * provides the implementation of main().
 */
//---------------------------------------------------------------------------//

#ifndef gtest_utils_gtest_hh
#define gtest_utils_gtest_hh

#include <Utils/config.h>
#include "gtest.h"
#include "Gtest_Functions.hh"

#include <iostream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"

// Import some common stuff into the global namespace
using std::cout;
using std::cerr;
using std::endl;
using profugus::soft_equiv;

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

//---------------------------------------------------------------------------//
// MAIN FUNCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Launch the google test harness from this unit test
 *
 * \warning We are intentionally putting a non-inlined function implementation
 * in the header! This is so that we can simply include "gtest/utils_gtest.hh"
 * and not have to rewrite this section of code.
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
#endif // gtest_utils_gtest_hh

//---------------------------------------------------------------------------//
//              end of gtest/utils_gtest.hh
//---------------------------------------------------------------------------//
