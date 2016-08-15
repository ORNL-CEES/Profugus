//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Atomic_Add_Kernel.cu
 * \author Seth R Johnson
 * \date   Thu Aug 15 11:09:56 2013
 * \brief  Atomic_Add_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../cuda_utils/Kernel_Header.cuh"

#include "Atomic_Add_Kernel.cuh"
#include "../cuda_utils/Atomic_Add.cuh"
#include "../cuda_utils/CudaDBC.hh"

namespace cuda_utils
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_T>
__global__ void atomic_add_test_kernel(
        Float_T* const __restrict__ out,
        unsigned int                num_increments)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (; i < num_increments; i += stride)
        Atomic_Add<Arch_Switch>::fetch_add(&out[0], static_cast<Float_T>(1));
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_Switch, typename Float_T>
void atomic_add_test(Atomic_Add_Kernel_Data<Arch_Switch, Float_T>& kd)
{
    REQUIRE(kd.output.size() == 1);
    REQUIRE(kd.num_increments > 0);

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
typedef ::cuda_utils::arch::Device Arch_t;
#else
typedef ::cuda_utils::arch::Host Arch_t;
#endif

template void atomic_add_test(Atomic_Add_Kernel_Data<Arch_t, float>& kd);
template void atomic_add_test(Atomic_Add_Kernel_Data<Arch_t, double>& kd);

//---------------------------------------------------------------------------//
} // end namespace cuda_utils

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Add_Kernel.cu
//---------------------------------------------------------------------------//
