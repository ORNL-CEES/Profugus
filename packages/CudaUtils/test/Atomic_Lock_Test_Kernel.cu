//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Atomic_Lock_Test_Kernel.cu
 * \author Seth R Johnson
 * \date   Thu Aug 15 08:16:56 2013
 * \brief  Atomic_Lock_Test_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../cuda_utils/Kernel_Header.cuh"

#include "Atomic_Lock_Test_Kernel.cuh"

#include "../cuda_utils/Atomic_Lock_Kernel.cuh"
#include "../cuda_utils/CudaDBC.hh"

namespace cuda_utils
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
__global__ void lock_test_kernel(
              Atomic_Lock_Kernel<Arch_Switch> lock,
              float* const __restrict__       out,
        unsigned int                          size)
{
    unsigned int i      = threadIdx.x;
    unsigned int stride = blockDim.x;

    // Acquire lock
    __syncthreads();
    if (threadIdx.x == 0)
        lock.acquire();
    __syncthreads();

    for (; i < size; i += stride)
        out[i] += 1;

    // Release lock
    __syncthreads();
    if (threadIdx.x == 0)
        lock.release();
    __syncthreads();
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
__host__ void lock_test(Lock_Kernel_Data<Arch_Switch>& kd)
{
    // Check for prior launch failure
    CudaCall(cudaGetLastError());

    CUDA_LAUNCH(lock_test_kernel<Arch_Switch>, kd.launch_args)(
            kd.lock.data(),
            kd.output.data(),
            kd.output.size());

    // Wait until kernel is finished
    CudaInsist(cudaDeviceSynchronize(), "Kernel execution error");
}

//---------------------------------------------------------------------------//
// INSTANTIATIONS
//---------------------------------------------------------------------------//
#ifdef __NVCC__
typedef ::cuda_utils::arch::Device Arch_t;
#else
typedef ::cuda_utils::arch::Host Arch_t;
#endif

template void lock_test(Lock_Kernel_Data<Arch_t>& kd);

//---------------------------------------------------------------------------//
} // end namespace cuda_utils

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Lock_Test_Kernel.cu
//---------------------------------------------------------------------------//
