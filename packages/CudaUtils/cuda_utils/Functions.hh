//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Functions.hh
 * \author Seth R Johnson
 * \date   Sat Sep 21 19:44:58 2013
 * \brief  Useful function declarations
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Functions_hh
#define CudaUtils_cuda_utils_Functions_hh

namespace cuda
{

// Return smallest value that satisfies result * divisor >= dividend
unsigned int ceil_quotient(unsigned int dividend, unsigned int divisor);

} // end namespace cuda

#endif // CudaUtils_cuda_utils_Functions_hh

//---------------------------------------------------------------------------//
//                 end of Functions.hh
//---------------------------------------------------------------------------//
