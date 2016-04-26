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
