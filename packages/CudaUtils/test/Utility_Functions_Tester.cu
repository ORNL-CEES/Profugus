//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Utility_Functions_Tester.cu
 * \author Tom Evans
 * \date   Mon Dec 05 14:42:30 2016
 * \brief  Utility_Functions_Tester member and kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>

#include "Utility_Functions_Tester.hh"

#include "cuda_utils/Utility_Functions.hh"

//---------------------------------------------------------------------------//
// THREAD ID TESTS
//---------------------------------------------------------------------------//

__global__
void kernel(int *tx,
            int *ty,
            int *tz,
            int *bx,
            int *by)
{
    auto tid = cuda::utility::thread_id();
    tx[tid]  = threadIdx.x;
    ty[tid]  = threadIdx.y;
    tz[tid]  = threadIdx.z;
    bx[tid]  = blockIdx.x;
    by[tid]  = blockIdx.y;
}

//---------------------------------------------------------------------------//

void ThreadID_Test::run(
    const unsigned int (&nb)[2],
    const unsigned int (&nt)[3])
{
    int num_threads = nt[0] * nt[1] * nt[2];
    int num_blocks  = nb[0] * nb[1];
    int N           = num_threads * num_blocks;

    thrust::device_vector<int> tx(N, -1), ty(N, -1), tz(N, -1);
    thrust::device_vector<int> bx(N, -1), by(N, -1);

    dim3 threads = {nt[0], nt[1], nt[2]};
    dim3 blocks  = {nb[0], nb[1]};

    kernel<<<blocks,threads>>>(tx.data().get(),
                               ty.data().get(),
                               tz.data().get(),
                               bx.data().get(),
                               by.data().get());

    // Make vectors on host and copy data into them from device
    vtx.resize(N);
    vty.resize(N);
    vtz.resize(N);
    vbx.resize(N);
    vby.resize(N);

    thrust::copy(tx.begin(), tx.end(), vtx.begin());
    thrust::copy(ty.begin(), ty.end(), vty.begin());
    thrust::copy(tz.begin(), tz.end(), vtz.begin());
    thrust::copy(bx.begin(), bx.end(), vbx.begin());
    thrust::copy(by.begin(), by.end(), vby.begin());
}

//---------------------------------------------------------------------------//

void ThreadID_Test::test_1D_1D()
{
    run({2,1}, {3,1,1});

    rtx = {0, 1, 2, 0, 1, 2};
    rty = {0, 0, 0, 0, 0, 0};
    rtz = {0, 0, 0, 0, 0, 0};

    rbx = {0, 0, 0, 1, 1, 1};
    rby = {0, 0, 0, 0, 0, 0};

    test();
}

//---------------------------------------------------------------------------//

void ThreadID_Test::test_1D_2D()
{
    run({2,1}, {3,2,1});

    rtx = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
    rty = {0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1};
    rtz = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    rbx = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
    rby = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    test();
}

//---------------------------------------------------------------------------//

void ThreadID_Test::test_1D_3D()
{
    run({2,1}, {3,2,2});

    rtx = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
           0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
    rty = {0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
           0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1};
    rtz = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
           0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};

    rbx = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    rby = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    test();
}

//---------------------------------------------------------------------------//

void ThreadID_Test::test_2D_1D()
{
    run({2,2}, {2,1,1});

    rtx = {0, 1, 0, 1, 0, 1, 0, 1};
    rty = {0, 0, 0, 0, 0, 0, 0, 0};
    rtz = {0, 0, 0, 0, 0, 0, 0, 0};

    rbx = {0, 0, 1, 1, 0, 0, 1, 1};
    rby = {0, 0, 0, 0, 1, 1, 1, 1};

    test();
}

//---------------------------------------------------------------------------//

void ThreadID_Test::test_2D_2D()
{
    run({2,2}, {2,2,1});

    rtx = {0, 1, 0, 1,
           0, 1, 0, 1,
           0, 1, 0, 1,
           0, 1, 0, 1};

    rty = {0, 0, 1, 1,
           0, 0, 1, 1,
           0, 0, 1, 1,
           0, 0, 1, 1};

    rtz = {0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0};

    rbx = {0, 0, 0, 0,
           1, 1, 1, 1,
           0, 0, 0, 0,
           1, 1, 1, 1};

    rby = {0, 0, 0, 0,
           0, 0, 0, 0,
           1, 1, 1, 1,
           1, 1, 1, 1};

    test();
}

//---------------------------------------------------------------------------//

void ThreadID_Test::test_2D_3D()
{
    run({2,2}, {2,2,2});

    rtx = {0, 1, 0, 1, 0, 1, 0, 1,
           0, 1, 0, 1, 0, 1, 0, 1,
           0, 1, 0, 1, 0, 1, 0, 1,
           0, 1, 0, 1, 0, 1, 0, 1};

    rty = {0, 0, 1, 1, 0, 0, 1, 1,
           0, 0, 1, 1, 0, 0, 1, 1,
           0, 0, 1, 1, 0, 0, 1, 1,
           0, 0, 1, 1, 0, 0, 1, 1};

    rtz = {0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1};

    rbx = {0, 0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 1, 1, 1, 1, 1,
           0, 0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 1, 1, 1, 1, 1};

    rby = {0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1};


    test();
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Utility_Functions_Tester.cu
//---------------------------------------------------------------------------//
