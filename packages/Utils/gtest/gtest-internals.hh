//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/gtest/gtest-internals.hh
 * \author Seth R Johnson
 * \date   Tue Apr 02 11:43:51 2013
 * \brief  Expose some google test internals that we want to use.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_gtest_gtest_internals_hh
#define Utils_gtest_gtest_internals_hh

//---------------------------------------------------------------------------//
// GTEST INTERNAL DECLARATIONS
//---------------------------------------------------------------------------//
// WARNING: this uses the google test internal API. If we update google test,
// make sure to update these declarations.
namespace testing {
namespace internal {

enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

void ColoredPrintf(GTestColor color, const char* fmt, ...);

} // end namespace internal
} // end namespace testing

//---------------------------------------------------------------------------//

#endif // Utils_gtest_gtest_internals_hh

//---------------------------------------------------------------------------//
//              end of gtest/gtest-internals.hh
//---------------------------------------------------------------------------//
