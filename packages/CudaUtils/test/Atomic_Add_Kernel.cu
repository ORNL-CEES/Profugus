//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Atomic_Add_Kernel.cu
 * \author Seth R Johnson
 * \date   Thu Aug 15 11:09:56 2013
 * \brief  Atomic_Add_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../Kernel_Header.cuh"

#include "Atomic_Add_Kernel.cuh"
#include "../Atomic_Add.cuh"
#include "../CudaDBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_T>
__global__ void atomic_add_test_kernel(
        Float_T* const __restrict__ out,
        unsigned int                num_increments)
{
    Atomic_Add<Arch_Switch, Float_T> atomic_add;

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (; i < num_increments; i += stride)
        atomic_add(&out[0], static_cast<Float_T>(1));
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename Float_T>
void atomic_add_test(Atomic_Add_Kernel_Data<Arch_Switch, Float_T>& kd)
{
    Require(kd.output.size() == 1);
    Require(kd.num_increments > 0);

    // Check for prior launch failure
    CudaCall(cudaGetLastError());

    CUDA_LAUNCH((atomic_add_test_kernel<Arch_Switch, Float_T>),
                kd.launch_args)(
            kd.output.data(),
            kd.num_increments);

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

template void atomic_add_test(Atomic_Add_Kernel_Data<Arch_t, float>& kd);
template void atomic_add_test(Atomic_Add_Kernel_Data<Arch_t, double>& kd);

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Add_Kernel.cu
//---------------------------------------------------------------------------//
