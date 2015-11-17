// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstPseudo_Cuda.cc
 * \author Seth R Johnson
 * \date   Tue Aug 13 14:09:29 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Pseudo_Cuda.hh"

//---------------------------------------------------------------------------//
// COMPILABILITY CHECKS
//---------------------------------------------------------------------------//

__device__ int some_pig()
{
    return 123;
}

__host__ void host_function()
{
    /* * */
}

__global__ void kernel_call(int* __restrict__ result)
{
    __shared__ float dummy[3];

    __syncthreads();
    *result = some_pig();
    __threadfence();
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(PseudoCudaTest, everything)
{
    int x;
    kernel_call(&x);
    EXPECT_EQ(123, x);

    // Check thread index, etc.
    EXPECT_EQ(1, gridDim.x);
    EXPECT_EQ(1, gridDim.y);
    EXPECT_EQ(1, gridDim.z);
    EXPECT_EQ(0, blockIdx.x);
    EXPECT_EQ(0, blockIdx.y);
    EXPECT_EQ(0, blockIdx.z);
    EXPECT_EQ(1, blockDim.x);
    EXPECT_EQ(1, blockDim.y);
    EXPECT_EQ(1, blockDim.z);
    EXPECT_EQ(0, threadIdx.x);
    EXPECT_EQ(0, threadIdx.y);
    EXPECT_EQ(0, threadIdx.z);

    EXPECT_EQ(1, warpSize);
}

//---------------------------------------------------------------------------//
//                        end of tstPseudo_Cuda.cc
//---------------------------------------------------------------------------//
