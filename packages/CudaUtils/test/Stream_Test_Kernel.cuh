//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Stream_Test_Kernel.cuh
 * \author Seth R Johnson
 * \date   Wed Oct 02 15:12:07 2013
 * \brief  Stream_Test_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Stream_Test_Kernel_cuh
#define cuda_utils_test_Stream_Test_Kernel_cuh

#include "Stream_Test_Kernel_Data.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_T, typename Float_T>
void stream_test(Stream_Test_Kernel_Data<Arch_T, Float_T>& kd);

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_test_Stream_Test_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Stream_Test_Kernel.cuh
//---------------------------------------------------------------------------//
