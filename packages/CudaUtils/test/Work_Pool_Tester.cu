//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Work_Pool_Tester.cu
 * \author Steven Hamilton
 * \date   Wed Mar 08 13:59:33 2017
 * \brief  Work_Pool_Tester kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/iterator/counting_iterator.h>
#include "Work_Pool_Tester.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Work_Pool.cuh"

//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

__global__ void iterate_kernel(cuda::Work_Pool pool)
{
    pool.initialize();
    
    int tid = cuda::utility::thread_id();
    int lane_id = tid % warpSize;
    int warp_id = tid / warpSize;

    if (lane_id == 0)
    {
        printf("Warp %i, original data\n",warp_id);
        printf("Warp %i, work per warp   = %i\n",warp_id,pool.work_per_warp());
        printf("Warp %i, work per thread = %i\n",warp_id,pool.work_per_thread());
    }
    __syncthreads();

    printf("On thread %i, first index = %i\n",tid,pool.work_id(0));
    printf("On thread %i, second index = %i\n",tid,pool.work_id(1));

    __syncthreads();

    // Inactivate some threads
    if (tid % 2 == 0)
        pool.set_inactive(0);
    if ((tid+1) % 2 == 0)
        pool.set_inactive(1);

    // Look at indices
    if (lane_id == 0)
    {
        printf("\n");
        printf("Warp %i, first modified data\n",warp_id);
        printf("Warp %i, work per warp   = %i\n",warp_id,pool.work_per_warp());
        printf("Warp %i, work per thread = %i\n",warp_id,pool.work_per_thread());
    }
    printf("On thread %i, first index = %i\n",tid,pool.work_id(0));
    printf("On thread %i, second index = %i\n",tid,pool.work_id(1));

    __syncthreads();

    // Consolidate
    pool.consolidate();

    // Get new indices
    if (lane_id == 0)
    {
        printf("\n");
        printf("Warp %i, first consolidated data\n",warp_id);
        printf("Warp %i, work per warp   = %i\n",warp_id,pool.work_per_warp());
        printf("Warp %i, work per thread = %i\n",warp_id,pool.work_per_thread());
    }
    printf("On thread %i, first index = %i\n",tid,pool.work_id(0));

    __syncthreads();

    // Disable some more threads
    if (lane_id >=8 && lane_id < 24)
        pool.set_inactive(0);

    // Get new indices
    if (lane_id == 0)
    {
        printf("\n");
        printf("Warp %i, second modified data\n",warp_id);
        printf("Warp %i, work per warp   = %i\n",warp_id,pool.work_per_warp());
        printf("Warp %i, work per thread = %i\n",warp_id,pool.work_per_thread());
    }

    __syncthreads();

    printf("On thread %i, first index = %i\n",tid,pool.work_id(0));

    __syncthreads();

    // Consolidate
    pool.consolidate();

    // Get new indices
    if (lane_id == 0)
    {
        printf("\n");
        printf("Warp %i, second consolidated data\n",warp_id);
        printf("Warp %i, work per warp   = %i\n",warp_id,pool.work_per_warp());
        printf("Warp %i, work per thread = %i\n",warp_id,pool.work_per_thread());
    }

    __syncthreads();

    printf("On thread %i, first index = %i\n",tid,pool.work_id(0));

}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void Work_Pool_Tester::test_pool()
{
    int num_warps = 1;
    int work_per_thread = 2;
    int total_work = num_warps * work_per_thread * 32;

    thrust::device_vector<int> indices;
    thrust::counting_iterator<int> it(0);
    indices.assign(it,it+total_work);

    std::cout << "Building work pool DMM" << std::endl;
    cuda::Work_Pool_DMM mgr(indices,num_warps);

    int block_size = std::min(num_warps*32,256);
    int num_blocks = total_work / block_size / work_per_thread;
    if (total_work % block_size)
        num_blocks++;
    std::cout << "Launching kernel: " << num_blocks << " " << block_size << std::endl;
    iterate_kernel<<<num_blocks,block_size>>>(mgr.device_instance());
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Work_Pool_Tester.cu
//---------------------------------------------------------------------------//
