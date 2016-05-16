//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Profiler_Kernel.cuh
 * \author Seth R Johnson
 * \date   Tue Jul 09 08:29:11 2013
 * \brief  Profiler_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Profiler_Kernel_cuh
#define cuda_utils_test_Profiler_Kernel_cuh

#include "../cuda_utils/Device_Vector.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_T, typename Float_T>
unsigned long int operation_test(
        unsigned int num_blocks,
        unsigned int num_threads,
        Device_Vector<Arch_T, Float_T>& data);

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_test_Profiler_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Profiler_Kernel.cuh
//---------------------------------------------------------------------------//
