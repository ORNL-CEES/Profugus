//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Polyglot_Kernel.cu
 * \author Seth R Johnson
 * \date   Tue Aug 13 15:32:33 2013
 * \brief  Polyglot_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../cuda_utils/Kernel_Header.cuh"

#include <utility>

#include "Polyglot_Kernel.cuh"
#include "../cuda_utils/CudaDBC.hh"
#include "../cuda_utils/Hardware.hh"
#include "../cuda_utils/Host_Vector.hh"

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
template<typename Arch_Switch>
__global__ void polyglot_copy_kernel_vector(
    const Device_Vector<Arch_Switch,float>* in,
    Device_Vector<Arch_Switch,float>* out )
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (; i < in->size(); i += stride)
        (*out)[i] = (*in)[i];
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
__host__ void polyglot_copy(Polyglot_Kernel_Data<Arch_Switch>& kd)
{
    REQUIRE(kd.input.size() == kd.output.size());

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
    REQUIRE(input.size() == output.size());

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
    REQUIRE(input.is_mapped());
    REQUIRE(input.size() == output.size());

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
void polyglot_copy_vector(
        const Device_Vector<Arch_Switch, Float_Type>& input,
              Device_Vector<Arch_Switch, Float_Type>& output)
{
    REQUIRE(input.size() == output.size());

    typedef Device_Vector<Arch_Switch,Float_Type> VecType;

    VecType* input_device;
    VecType* output_device;

#ifdef __NVCC__
    cudaMalloc( (void**) &input_device, sizeof(VecType) );
    cudaMemcpy( input_device, &input, sizeof(VecType),
		cudaMemcpyHostToDevice );

    cudaMalloc( (void**) &output_device, sizeof(VecType) );
    cudaMemcpy( output_device, &output, sizeof(VecType),
		cudaMemcpyHostToDevice );

#else
    input_device = const_cast<VecType*>(&input);
    output_device = &output;
#endif

    typedef ::cuda::Hardware<Arch_Switch> Hardware_t;

    unsigned int num_threads = Hardware_t::num_cores_per_mp();
    unsigned int num_blocks  = Hardware_t::num_multiprocessors();

    polyglot_copy_kernel_vector
#ifdef __NVCC__
        <<<num_blocks, num_threads>>>
#endif
        ( input_device, output_device );

    // Wait until kernel is finished
    CudaInsist(cudaDeviceSynchronize(), "Kernel execution error");

#ifdef __NVCC__
    cudaFree( input_device );
    cudaFree( output_device );    
#endif
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

template void polyglot_copy_vector(
        const Device_Vector<Arch_t, float>&, Device_Vector<Arch_t, float>&);


//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Polyglot_Kernel.cu
//---------------------------------------------------------------------------//
