//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Functions.cc
 * \author Seth R Johnson
 * \date   Sun Sep 22 11:08:49 2013
 * \brief  Useful function definitions
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Functions.hh"

#include "harness/DBC.hh"

namespace cuda_utils
{
//---------------------------------------------------------------------------//
/*!
 * \brief Round up to the nearest divisor, then divide
 *
 * This is useful for calculating block counts that, given a number of
 * threads per block (divisor) and a total thread count (dividend), ensures
 * that result * divisor >= dividend
 *
 * \param dividend Integer to be divided
 * \param divisor
 */
unsigned int ceil_quotient(unsigned int dividend, unsigned int divisor)
{
    REQUIRE(divisor > 0);
    unsigned int result = (dividend + divisor - 1) / divisor;
    ENSURE(result * divisor >= dividend);
    return result;
}

//---------------------------------------------------------------------------//

} // end namespace cuda_utils

//---------------------------------------------------------------------------//
//                 end of Functions.cc
//---------------------------------------------------------------------------//
