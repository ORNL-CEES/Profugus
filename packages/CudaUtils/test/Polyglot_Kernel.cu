//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Polyglot_Kernel.cu
 * \author Seth R Johnson
 * \date   Tue Aug 13 15:32:33 2013
 * \brief  Polyglot_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../Kernel_Header.cuh"

#include <utility>

#include "Polyglot_Kernel.cuh"
#include "../CudaDBC.hh"
#include "../Hardware.hh"
#include "../Host_Vector.hh"

//---------------------------------------------------------------------------//
// DEFINITIONS
//---------------------------------------------------------------------------//
namespace cuda
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
__global__ void polyglot_copy_kernel(
        const float* const __restrict__ in,
              float* const __restrict__ out,
        unsigned int                    size)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (; i < size; i += stride)
        out[i] = in[i];
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
__host__ void polyglot_copy(Polyglot_Kernel_Data<Arch_Switch>& kd)
{
    Require(kd.input.size() == kd.output.size());

    // Check for prior launch failure
    CudaCall(cudaGetLastError());

    CUDA_LAUNCH(polyglot_copy_kernel<Arch_Switch>, kd.launch_args)(
            kd.input.data(),
            kd.output.data(),
            kd.output.size());

    // Wait until kernel is finished
    CudaInsist(cudaDeviceSynchronize(), "Kernel execution error");
}

//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_Type>
void polyglot_copy(
        const Device_Vector<Arch_Switch, Float_Type>& input,
              Device_Vector<Arch_Switch, Float_Type>& output)
{
    Require(input.size() == output.size());

    typedef ::cuda::Hardware<Arch_Switch> Hardware_t;

    unsigned int num_threads = Hardware_t::num_cores_per_mp();
    unsigned int num_blocks  = Hardware_t::num_multiprocessors();

    polyglot_copy_kernel<Arch_Switch>
#ifdef __NVCC__
        <<<num_blocks, num_threads>>>
#endif
        (
            input.data(),
            output.data(),
            output.size());

    // Wait until kernel is finished
    CudaInsist(cudaDeviceSynchronize(), "Kernel execution error");
}

//---------------------------------------------------------------------------//
template<typename Arch_Switch, typename Float_Type>
void polyglot_copy(
        const Host_Vector<Float_Type>&   input,
              Device_Vector<Arch_Switch,Float_Type>& output)
{
    Require(input.is_mapped());
    Require(input.size() == output.size());

    typedef ::cuda::Hardware<Arch_Switch> Hardware_t;

    unsigned int num_threads = Hardware_t::num_cores_per_mp();
    unsigned int num_blocks  = Hardware_t::num_multiprocessors();

    polyglot_copy_kernel<Arch_Switch>
#ifdef __NVCC__
        <<<num_blocks, num_threads>>>
#endif
        (
            input.data(),
            output.data(),
            output.size());

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

template void polyglot_copy(Polyglot_Kernel_Data<Arch_t>& kd);

template void polyglot_copy(
        const Device_Vector<Arch_t, float>&, Device_Vector<Arch_t, float>&);

template void polyglot_copy(
        const Host_Vector<float>&, Device_Vector<Arch_t, float>&);

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Polyglot_Kernel.cu
//---------------------------------------------------------------------------//
