//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/gtest/Gtest_Functions.hh
 * \author Seth R Johnson
 * \date   Tue Apr 02 11:41:41 2013
 * \brief  Gtest_Functions class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_gtest_Gtest_Functions_hh
#define Utils_gtest_Gtest_Functions_hh

#include <cstddef>

namespace profugus
{
//===========================================================================//

// Implementation for the harness' "main" function.
int gtest_main(int argc, char *argv[]);

// Print the "skip" message from the skip macro
void print_skip_message();

// Calculate the number of digits in a large number
unsigned int calc_num_digits(std::size_t number);

//===========================================================================//
} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Gtest_Functions.i.hh"
//---------------------------------------------------------------------------//
#endif // Utils_gtest_Gtest_Functions_hh

//---------------------------------------------------------------------------//
//              end of gtest/Gtest_Functions.hh
//---------------------------------------------------------------------------//
