//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Atomic_Add_Kernel.cuh
 * \author Seth R Johnson
 * \date   Thu Aug 15 11:09:56 2013
 * \brief  Atomic_Add_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Atomic_Add_Kernel_cuh
#define cuda_utils_test_Atomic_Add_Kernel_cuh

#include "Atomic_Add_Kernel_Data.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename Float_T>
void atomic_add_test(Atomic_Add_Kernel_Data<Arch_Switch, Float_T>& kd);

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_test_Atomic_Add_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Add_Kernel.cuh
//---------------------------------------------------------------------------//
