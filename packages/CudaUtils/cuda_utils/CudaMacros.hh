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

#include "CudaDBC.hh"

//---------------------------------------------------------------------------//
// Prefix for device-only functions.
#ifdef __NVCC__
#define PROFUGUS_DEVICE_FUNCTION __device__ 
#else
#define PROFUGUS_DEVICE_FUNCTION
#endif

// Prefix for host-device functions.
#ifdef __NVCC__
#define PROFUGUS_HOST_DEVICE_FUNCTION __host__ __device__
#else
#define PROFUGUS_HOST_DEVICE_FUNCTION
#endif

// Insist only on the device.
#ifdef __NVCC__
#define PROFUGUS_INSIST_ON_DEVICE
#else
#define PROFUGUS_INSIST_ON_DEVICE \
    do { INSIST(false,"Calling device code from host!"); } while(0)
#endif

//---------------------------------------------------------------------------//

#endif // cuda_utils_CudaMacros_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/CudaMacros.hh
//---------------------------------------------------------------------------//
