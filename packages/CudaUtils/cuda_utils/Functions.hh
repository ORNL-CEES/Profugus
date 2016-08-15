//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Functions.hh
 * \author Seth R Johnson
 * \date   Sat Sep 21 19:44:58 2013
 * \brief  Useful function declarations
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Functions_hh
#define cuda_utils_Functions_hh

namespace cuda_utils
{

// Return smallest value that satisfies result * divisor >= dividend
unsigned int ceil_quotient(unsigned int dividend, unsigned int divisor);

} // end namespace cuda_utils

#endif // cuda_utils_Functions_hh

//---------------------------------------------------------------------------//
//                 end of Functions.hh
//---------------------------------------------------------------------------//
