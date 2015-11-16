//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Stream_Test_Kernel.cu
 * \author Seth R Johnson
 * \date   Wed Oct 02 15:12:07 2013
 * \brief  Stream_Test_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../Kernel_Header.cuh"

#include "Stream_Test_Kernel.cuh"

#include "../CudaDBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

template<typename Arch_T, typename float_type>
__global__ void stream_test_kernel(
        const float_type* const __restrict__ all_tau,
        const float_type* const __restrict__ all_src,
        const float_type* const __restrict__ input,
              float_type* const __restrict__ output,
        const unsigned int                   num_rays,
        const unsigned int                   num_segments)
{
    unsigned int ray    = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (; ray < num_rays; ray += stride)
    {
        // Get input value
        float_type psi = input[ray];

        // Pointer for traversing thicknesses
        const float_type* __restrict__ tau_it = all_tau + ray;
        const float_type* __restrict__ src_it = all_src + ray;

        for (unsigned int seg = 0; seg < num_segments;
                ++seg, tau_it += num_rays, src_it += num_rays)
        {
            Check(tau_it < all_tau + num_rays * num_segments);
            Check(src_it < all_src + num_rays * num_segments);

            const float_type tau = *tau_it;
            const float_type src = *src_it;

            Check(tau >= 0 && tau < 1);
            Check(src >= 0);

            // Do the segments multiple times for fun
            for (unsigned int i = 0; i < 1024; ++i)
            {
                // Approximate taylor whatsit for (1 - exp(-tau)) for tau ~ 0
                psi = src + psi * (tau
                                    * (1 - tau / 2
                                            * (1 - tau / 3
                                                    * (1 - tau / 4))));
            }
        }

        // Save output value
        output[ray] = psi;
    }
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

template<typename Arch_T, typename Float_T>
void stream_test(Stream_Test_Kernel_Data<Arch_T, Float_T>& kd)
{
    Require(kd.input.size() == kd.output.size());
    Require(kd.input.size() == kd.num_rays);
    Require(kd.tau.size() == kd.num_rays * kd.num_segments);

    CUDA_LAUNCH((stream_test_kernel<Arch_T, Float_T>), kd.launch_args)(
            kd.tau.cdata(),
            kd.source.cdata(),
            kd.input.cdata(),
            kd.output.data(),
            kd.num_rays,
            kd.num_segments);
}

//---------------------------------------------------------------------------//
// INSTANTIATIONS
//---------------------------------------------------------------------------//
#ifdef __NVCC__
typedef ::cuda::arch::Device Arch_t;
#else
typedef ::cuda::arch::Host Arch_t;
#endif

template void stream_test(Stream_Test_Kernel_Data<Arch_t,float>&);
template void stream_test(Stream_Test_Kernel_Data<Arch_t,double>&);


//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Stream_Test_Kernel.cu
//---------------------------------------------------------------------------//
