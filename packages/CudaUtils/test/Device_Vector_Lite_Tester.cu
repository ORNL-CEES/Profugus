//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/Device_Vector_Lite_Tester.cu
 * \author Steven Hamilton
 * \date   Tue May 03 12:09:58 2016
 * \brief  Device_Vector_Lite_Tester class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "Device_Vector_Lite_Tester.hh"

#include "gtest/Gtest_Functions.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Device_Vector_Lite.hh"

using Vec_Lite      = cuda_utils::Device_Vector_Lite<int, 3>;
using Host_Vec_Lite = Vec_Lite::Host_Vec_Lite;

//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

__global__ void from_device(Vec_Lite *vals, int num_vals)
{
    DEVICE_CHECK(num_vals == 2);

    int tid = threadIdx.x + blockIdx.x*blockDim.x;

    if (tid == 0)
        vals[tid] = {2, 4, 0};
    else if (tid == 1)
        vals[tid] = {1, 3, 5};
}

//---------------------------------------------------------------------------//

__global__ void copy_back(const Vec_Lite *vals_in,
                                Vec_Lite *vals_out,
                                int       num_vals)
{
    DEVICE_CHECK(num_vals == 2);

    int tid = threadIdx.x + blockIdx.x*blockDim.x;
    if (tid < num_vals)
    {
        Vec_Lite tmp;
        for (auto i : {0, 1, 2})
            tmp[i] = vals_in[tid][i];
        vals_out[tid] = tmp;
    }
}

//---------------------------------------------------------------------------//
// FUNCTIONS
//---------------------------------------------------------------------------//

void Device_Vector_Lite_Tester::test_host()
{
    // Test host API

    // Check adherence to concepts
    EXPECT_TRUE(std::is_pod<Vec_Lite>::value);
    EXPECT_TRUE(std::is_standard_layout<Vec_Lite>::value);

    Vec_Lite a = {1, 4, -3};
    EXPECT_EQ(1, a[0]);
    EXPECT_EQ(4, a[1]);
    EXPECT_EQ(-3,a[2]);

    Vec_Lite b;
    b[0] = 2;
    b[1] = 0;
    b[2] = 12;
    EXPECT_EQ(2, b[0]);
    EXPECT_EQ(0, b[1]);
    EXPECT_EQ(12,b[2]);
}

//---------------------------------------------------------------------------//

void Device_Vector_Lite_Tester::test_device()
{
    // Test copying from device
    int num_vals = 2;
    thrust::device_vector<Vec_Lite> v1(num_vals);

    from_device<<<1,num_vals>>>(v1.data().get(), num_vals);

    thrust::host_vector<Vec_Lite> host_v1 = v1;
    EXPECT_EQ(2, host_v1[0][0]);
    EXPECT_EQ(4, host_v1[0][1]);
    EXPECT_EQ(0, host_v1[0][2]);
    EXPECT_EQ(1, host_v1[1][0]);
    EXPECT_EQ(3, host_v1[1][1]);
    EXPECT_EQ(5, host_v1[1][2]);

    // Make sure we can perform copy from one array to another on device
    thrust::device_vector<Vec_Lite> v2(num_vals);
    copy_back<<<1,num_vals>>>(v1.data().get(),v2.data().get(),num_vals);

    thrust::host_vector<Vec_Lite> host_v2 = v2;
    EXPECT_EQ(2, host_v2[0][0]);
    EXPECT_EQ(4, host_v2[0][1]);
    EXPECT_EQ(0, host_v2[0][2]);
    EXPECT_EQ(1, host_v2[1][0]);
    EXPECT_EQ(3, host_v2[1][1]);
    EXPECT_EQ(5, host_v2[1][2]);

    // And a straight on-device copy using thrust
    thrust::device_vector<Vec_Lite> v3 = v1;
    thrust::host_vector<Vec_Lite> host_v3 = v3;
    EXPECT_EQ(2, host_v3[0][0]);
    EXPECT_EQ(4, host_v3[0][1]);
    EXPECT_EQ(0, host_v3[0][2]);
    EXPECT_EQ(1, host_v3[1][0]);
    EXPECT_EQ(3, host_v3[1][1]);
    EXPECT_EQ(5, host_v3[1][2]);
}

//---------------------------------------------------------------------------//

void Device_Vector_Lite_Tester::test_host_copy()
{
    Host_Vec_Lite h = {1, 2, 3};
    Vec_Lite      d = {0, 0, 0};

    EXPECT_EQ(1, h[0]);
    EXPECT_EQ(2, h[1]);
    EXPECT_EQ(3, h[2]);

    EXPECT_EQ(0, d[0]);
    EXPECT_EQ(0, d[1]);
    EXPECT_EQ(0, d[2]);

    d = Vec_Lite::from_host(h);

    h[2] = 5;

    EXPECT_EQ(1, h[0]);
    EXPECT_EQ(2, h[1]);
    EXPECT_EQ(5, h[2]);

    EXPECT_EQ(1, d[0]);
    EXPECT_EQ(2, d[1]);
    EXPECT_EQ(3, d[2]);
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Device_Vector_Lite_Tester.cu
//---------------------------------------------------------------------------//
