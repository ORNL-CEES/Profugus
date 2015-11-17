//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/CudaDBC.hh
 * \author Seth R Johnson
 * \date   Thu Jun 27 15:21:33 2013
 * \brief  DBC class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_CudaDBC_hh
#define cuda_utils_CudaDBC_hh

#include "config.h"

//---------------------------------------------------------------------------//
// If compiling device code, disable DBC
//---------------------------------------------------------------------------//
#ifdef __CUDA_ARCH__
#ifdef Require
// Undefine DBC macros since we've included DBC.hh
#include "harness/DBC_undef.hh"
#else // Require
// Prevent DBC from being loaded
#define harness_DBC_hh
#endif // Require

// Null-ops for DBC macros
#include "harness/DBC_nulldef.hh"
#else // __CUDA_ARCH__
#include "harness/DBC.hh"
#endif

//---------------------------------------------------------------------------//
// Compile assertions on host code when CUDA is enabled
#if defined(USE_CUDA) && !defined(__CUDA_ARCH__) && !defined(PSEUDO_CUDA)

//---------------------------------------------------------------------------//
// FUNCTION DEFINITIONS
//---------------------------------------------------------------------------//
namespace cuda
{

// Raise profugus::assertion with a formatted error message
void toss_cuda_cookies(
        const char* err_string,
        const char* msg,
        const char* file,
        int         line);

} // end namespace cuda

//---------------------------------------------------------------------------//
// MACRO DEFINITIONS
//---------------------------------------------------------------------------//

/*!
 * \def CudaCall(statement)
 *
 * Execute the wrapped statement and throw a message if it fails.
 *
 * If it fails, we call \c cudaGetLastError to clear the error code.
 *
 * \code
 *     CudaCall(cudaMalloc(&ptr_gpu, 100 * sizeof(float)));
 * \endcode
 */
#define CudaCall(COND) \
    do \
    { \
        cudaError_t result_ = (COND); \
        if (result_ != cudaSuccess) \
        { \
            cudaGetLastError(); \
            ::cuda::toss_cuda_cookies(cudaGetErrorString(result_), #COND, \
                __FILE__, __LINE__); \
        } \
    } while (0);

/*!
 * \def CudaInsist(statement/condition, message)
 *
 * Call the given function and throw a message if it doesn't return cudaSuccess.
 *
 * If it fails, we call \c cudaGetLastError to clear the error code.
 */
#define CudaInsist(COND, MSG) \
    do \
    { \
        cudaError_t result_ = (COND); \
        if (result_ != cudaSuccess) \
        { \
            cudaGetLastError(); \
            ::cuda::toss_cuda_cookies(MSG, cudaGetErrorString(result_), \
                __FILE__, __LINE__); \
        } \
    } while (0);

//---------------------------------------------------------------------------//
#else // CUDA is disabled or we're compiling device code

// Compile out the CUDA calls
#define CudaInsist(COND, MSG)
#define CudaCall(COND)

#endif
//---------------------------------------------------------------------------//

#endif // cuda_utils_CudaDBC_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/CudaDBC.hh
//---------------------------------------------------------------------------//
