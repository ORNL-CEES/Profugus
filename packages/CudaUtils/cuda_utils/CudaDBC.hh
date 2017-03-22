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
#include "Utils/config.h"

//
// On-Device DBC using "assert"
//

#ifndef __APPLE__

#include <assert.h>

/* The following definitions of assertions help ensure that they're always
 * followed by a semicolon, and that if assertions are disabled we don't get
 * "unused variable" warnings.
 *
 * See:
 * http://cnicholson.net/2009/02/stupid-c-tricks-adventures-in-assert/
 */
#define UTILS_NOASSERT_(COND) \
    do { (void)sizeof(COND); } while (0)

// Insist is always on
// If condition fails, print message then assert
// Don't expect fancy stream machinery to work
//  (e.g. INSIST(n>0,"Value is " << n); ) like it would in normal DBC
#define DEVICE_INSIST(COND,MSG) \
    do \
    { \
        if (!(COND)) \
        { \
            printf(MSG "\n");  \
        } \
        assert(COND); \
    } while (0)

#if UTILS_DBC & 1
#define DEVICE_REQUIRE(COND) \
    do { assert(COND); } while (0)
#else
#define DEVICE_REQUIRE(COND) UTILS_NOASSERT_(COND)
#endif

#if UTILS_DBC & 2
#define DEVICE_CHECK(COND) \
    do { assert(COND); } while (0)
#else
#define DEVICE_CHECK(COND) UTILS_NOASSERT_(COND)
#endif

#if UTILS_DBC & 4
#define DEVICE_ENSURE(COND) \
    do { assert(COND); } while (0)
#else
#define DEVICE_ENSURE(COND) UTILS_NOASSERT_(COND)
#endif

#else   // __APPLE__
// No device-side asserts on Mac, good luck
#define DEVICE_INSIST(COND,MSG) UTILS_NOASSERT_(COND)
#define DEVICE_REQUIRE(COND)    UTILS_NOASSERT_(COND)
#define DEVICE_CHECK(COND)      UTILS_NOASSERT_(COND)
#define DEVICE_ENSURE(COND)     UTILS_NOASSERT_(COND)
#endif  // __APPLE__

//---------------------------------------------------------------------------//
// Compile assertions on host code when CUDA is enabled
#if defined(USE_CUDA) && !defined(__CUDA_ARCH__) && !defined(PSEUDO_CUDA)

//---------------------------------------------------------------------------//
// FUNCTION DEFINITIONS
//---------------------------------------------------------------------------//
namespace cuda_utils
{

// Raise profugus::assertion with a formatted error message
void toss_cuda_cookies(
        const char* err_string,
        const char* msg,
        const char* file,
        int         line);

} // end namespace cuda_utils

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
            ::cuda_utils::toss_cuda_cookies(cudaGetErrorString(result_), #COND, \
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
            ::cuda_utils::toss_cuda_cookies(MSG, cudaGetErrorString(result_), \
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
