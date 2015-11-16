//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Exp_Kernel.cu
 * \author Seth R Johnson
 * \date   Fri Aug 16 10:16:13 2013
 * \brief  Exp_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../Kernel_Header.cuh"

#include "Exp_Kernel.cuh"

#include "../CudaDBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_T>
__global__ void exp_test_kernel(
        Float_T* const __restrict__ out,
        unsigned int                size,
        Exponential<Arch_Switch, Float_T> exp)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (; i < size; i += stride)
        out[i] = exp(out[i]);
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename Float_T>
void exp_test(Device_Vector<Arch_Switch, Float_T>& inout,
        Exponential<Arch_Switch, Float_T> exp)
{
    Require(inout.size() == 1024);

    // Check for prior launch failure
    CudaCall(cudaGetLastError());

    exp_test_kernel<Arch_Switch, Float_T>
#ifdef __NVCC__
        <<<32, 32>>>
#endif
        (
            inout.data(),
            inout.size(),
            exp);

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

template void exp_test(Device_Vector<Arch_t,float>& ,
        Exponential<Arch_t,float> );
template void exp_test(Device_Vector<Arch_t,double>&,
        Exponential<Arch_t,double>);

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Exp_Kernel.cu
//---------------------------------------------------------------------------//
