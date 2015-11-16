//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Texture_Vector_Test_Kernel.cuh
 * \author Seth R Johnson
 * \date   Sat Sep 21 11:59:41 2013
 * \brief  Texture_Vector_Test_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Texture_Vector_Test_Kernel_cuh
#define cuda_utils_test_Texture_Vector_Test_Kernel_cuh

#include "../Texture_Vector.hh"
#include "../Device_Vector.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename T>
void texture_vector_test(
        const Texture_Vector<Arch_Switch, T>& input,
              Device_Vector<Arch_Switch, T>&  output);

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_test_Texture_Vector_Test_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Texture_Vector_Test_Kernel.cuh
//---------------------------------------------------------------------------//
