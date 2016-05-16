//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Kernel_Header.cuh
 * \author Seth R Johnson
 * \date   Wed Nov 27 09:22:28 2013
 * \brief  Kernel_Header kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Kernel_Header_cuh
#define cuda_utils_Kernel_Header_cuh

//---------------------------------------------------------------------------//
// Include this file at the TOP of all .cu files to reduce compile time.
//---------------------------------------------------------------------------//

#ifdef __CUDA_ARCH__

// Make sure DBC doesn't get included (don't compile in iostream, etc.)
#define harness_DBC_hh
#include "harness/DBC_nulldef.hh"

#endif // __CUDA_ARCH__

//---------------------------------------------------------------------------//

#endif // cuda_utils_Kernel_Header_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Kernel_Header.cuh
//---------------------------------------------------------------------------//
