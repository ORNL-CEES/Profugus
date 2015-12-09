//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Launch_Args.t.cuh
 * \author Seth R Johnson
 * \date   Wed Oct 02 13:16:37 2013
 * \brief  Launch class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Launch_Args_t_cuh
#define cuda_utils_Launch_Args_t_cuh

#include "Launch_Args.hh"
#include "cuda_runtime.h"

namespace cuda
{
//---------------------------------------------------------------------------//
// GLOBAL CUDA KERNEL
//---------------------------------------------------------------------------//
template<class Kernel>
__global__ void cuda_kernel( Kernel kernel, std::size_t N )
{
    std::size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if( idx < N )
        kernel( idx );
}

//---------------------------------------------------------------------------//
// PARALLEL LAUNCH
//---------------------------------------------------------------------------//
// Cuda specialization.
template <class Kernel>
void parallel_launch(
    Kernel& kernel, const Launch_Args<cuda::arch::Device>& launch_args )
{
    REQUIRE( launch_args.is_valid() );
    cuda_kernel<<<launch_args.grid_size(),
	launch_args.block_size(),
	launch_args.shared_mem(),
	launch_args.stream_handle()>>>( kernel, launch_args.num_elements() );
}

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_Launch_Args_t_cuh

//---------------------------------------------------------------------------//
//                 end of Launch_Args.t.cuh
//---------------------------------------------------------------------------//
