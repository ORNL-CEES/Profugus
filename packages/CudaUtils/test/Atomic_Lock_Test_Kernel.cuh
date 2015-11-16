//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Atomic_Lock_Test_Kernel.cuh
 * \author Seth R Johnson
 * \date   Thu Aug 15 08:16:56 2013
 * \brief  Atomic_Lock_Test_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Lock_Kernel_cuh
#define cuda_utils_test_Lock_Kernel_cuh

#include "Lock_Kernel_Data.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch>
void lock_test(Lock_Kernel_Data<Arch_Switch>& kd);

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_test_Lock_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Lock_Test_Kernel.cuh
//---------------------------------------------------------------------------//
