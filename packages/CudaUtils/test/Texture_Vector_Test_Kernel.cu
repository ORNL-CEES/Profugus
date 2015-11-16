//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Texture_Vector_Test_Kernel.cu
 * \author Seth R Johnson
 * \date   Sat Sep 21 11:59:41 2013
 * \brief  Texture_Vector_Test_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../Kernel_Header.cuh"

#include "Texture_Vector_Test_Kernel.cuh"

#include "../CudaDBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename T>
__global__ void texture_vector_test_kernel(
        const Texture_Vector_Kernel<Arch_Switch, T> input,
        size_t size,
        T*  output)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    for (; i < size; i += blockDim.x * gridDim.x)
    {
        output[i] = input[i];
    }
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename T>
void texture_vector_test(
        const Texture_Vector<Arch_Switch, T>& input,
              Device_Vector<Arch_Switch, T>&  output)
{
    Require(input.size() == output.size());

    // Check for prior launch failure
    CudaCall(cudaGetLastError());

    unsigned int num_threads = 64;
    unsigned int num_blocks  = output.size() / 64;

    texture_vector_test_kernel<Arch_Switch, T>
#ifdef __NVCC__
        <<<num_blocks, num_threads>>>
#endif
        (
            input.data(),
            input.size(),
            output.data());

    // Wait until kernel is finished
    CudaInsist(cudaDeviceSynchronize(), "Kernel execution error");
}

//---------------------------------------------------------------------------//
// INSTANTIATIONS
//---------------------------------------------------------------------------//
#ifdef __NVCC__
typedef ::cuda::arch::Device Arch_t;
#else
typedef ::cuda::arch::Host Arch_t;
#endif

template void texture_vector_test(
        const Texture_Vector<Arch_t, int>&,
              Device_Vector<Arch_t, int>&);
template void texture_vector_test(
        const Texture_Vector<Arch_t, float>&,
              Device_Vector<Arch_t, float>&);
template void texture_vector_test(
        const Texture_Vector<Arch_t, double>&,
              Device_Vector<Arch_t, double>&);


//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Texture_Vector_Test_Kernel.cu
//---------------------------------------------------------------------------//
