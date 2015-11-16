//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Exp_Kernel.cuh
 * \author Seth R Johnson
 * \date   Fri Aug 16 10:16:13 2013
 * \brief  Exp_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Exp_Kernel_cuh
#define cuda_utils_test_Exp_Kernel_cuh

#include "../Device_Vector.hh"
#include "../Exponential.cuh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename Float_T>
void exp_test(
        Device_Vector<Arch_Switch, Float_T>& inout,
        Exponential<Arch_Switch, Float_T>    exp);

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_test_Exp_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Exp_Kernel.cuh
//---------------------------------------------------------------------------//
