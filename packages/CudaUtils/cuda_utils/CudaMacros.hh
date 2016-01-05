//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/CudaMacros.hh
 * \author Stuart Slattery
 * \brief  Function macros for cuda programming.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_CudaMacros_hh
#define cuda_utils_CudaMacros_hh

//---------------------------------------------------------------------------//
// Wrapper for device-only functions. Will only be visible to NVCC compiler.
#ifdef __NVCC__
#define PROFUGUS_DEVICE_FUNCTION(...) __device__ __VA_ARGS__
#else
#define PROFUGUS_DEVICE_FUNCTION(...)
#endif

// Wrapper for host-device functions. Will be visible to both host compiler
// and NVCC compiler.
#ifdef __NVCC__
#define PROFUGUS_HOST_DEVICE_FUNCTION(...) __host__ __device__ __VA_ARGS__
#else
#define PROFUGUS_HOST_DEVICE_FUNCTION(...) __VA_ARGS__
#endif

//---------------------------------------------------------------------------//

#endif // cuda_utils_CudaMacros_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/CudaMacros.hh
//---------------------------------------------------------------------------//
