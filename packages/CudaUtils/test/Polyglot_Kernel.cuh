//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Polyglot_Kernel.cuh
 * \author Seth R Johnson
 * \date   Tue Aug 13 15:32:33 2013
 * \brief  Polyglot_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Polyglot_Kernel_cuh
#define cuda_utils_test_Polyglot_Kernel_cuh

#include "Polyglot_Kernel_Data.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
void polyglot_copy(Polyglot_Kernel_Data<Arch_Switch>& data);

//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_Type>
void polyglot_copy(
        const Device_Vector<Arch_Switch, Float_Type>& input,
              Device_Vector<Arch_Switch, Float_Type>& output);

//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_Type>
void polyglot_copy(
        const cuda::Host_Vector<Float_Type>&   in,  // Must be mapped memory
              cuda::Device_Vector<Arch_Switch,Float_Type>& out);

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // cuda_utils_test_Polyglot_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Polyglot_Kernel.cuh
//---------------------------------------------------------------------------//
