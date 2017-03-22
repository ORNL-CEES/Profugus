//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Device_View_Field_Tester.cu
 * \author Steven Hamilton
 * \date   Fri Nov 18 11:55:55 2016
 * \brief  Device_View_Field_Tester kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>

#include "gtest/Gtest_Functions.hh"
#include "cuda_utils/Device_View_Field.hh"
#include "Device_View_Field_Tester.hh"

typedef cuda::Device_View_Field<int>       DVF;
typedef cuda::const_Device_View_Field<int> const_DVF;

//---------------------------------------------------------------------------//
// GLOBAL FUNCTIONS
//---------------------------------------------------------------------------//

__global__ void write_to_view(const int *vals, DVF view)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    DEVICE_REQUIRE(tid < view.size());

    view[tid] = vals[tid];
}

__global__ void write_from_view(DVF view, int *vals)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    DEVICE_REQUIRE(tid < view.size());

    const_DVF new_view = view;

    vals[tid] = new_view[tid];
}

__global__ void write_from_const_view(const_DVF view, int *vals)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    DEVICE_REQUIRE(tid < view.size());

    const_DVF new_view = view;

    vals[tid] = new_view[tid];
}

__global__ void check_null(DVF view, int *vals)
{
    DEVICE_REQUIRE(vals[0] == -1);
    DEVICE_REQUIRE(vals[1] == -1);

    vals[0] = static_cast<int>(view.empty());
    vals[1] = view.size();
}

//---------------------------------------------------------------------------//
// TEST FUNCTIONS
//---------------------------------------------------------------------------//

void Device_View_Field_Tester::test_views()
{
    std::vector<int> ref = {6, 4, 2, 0, 1, 3, 5, 7};
    thrust::device_vector<int> int_vec = ref;

    thrust::device_vector<int> result(8,-1);

    // Write data from raw pointer into View
    write_to_view<<<1,8>>>(int_vec.data().get(),cuda::make_view(result));

    std::vector<int> host_result(8,-1);
    thrust::copy(result.begin(),result.end(),host_result.begin());
    EXPECT_VEC_EQ(ref,host_result);

    // Write data from view into raw pointer
    result.resize(8,-1);
    write_from_view<<<1,8>>>(cuda::make_view(int_vec),
                             result.data().get());

    host_result.resize(8,-1);
    thrust::copy(result.begin(),result.end(),host_result.begin());
    EXPECT_VEC_EQ(ref,host_result);

    // Write data from const view into raw pointer
    const thrust::device_vector<int> const_int_vec = int_vec;
    result.resize(8,-1);
    write_from_const_view<<<1,8>>>(cuda::make_view(const_int_vec),
                                   result.data().get());

    host_result.resize(8,-1);
    thrust::copy(result.begin(),result.end(),host_result.begin());
    EXPECT_VEC_EQ(ref,host_result);

    // Make an empty device vector
    thrust::device_vector<int> null;
    EXPECT_TRUE(null.empty());

    result.resize(2);
    thrust::fill(result.begin(), result.end(), -1);
    check_null<<<1,1>>>(cuda::make_view(null), result.data().get());

    host_result.resize(2,-1);
    thrust::copy(result.begin(),result.end(),host_result.begin());

    EXPECT_EQ(1, host_result[0]);
    EXPECT_EQ(0, host_result[1]);
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Device_View_Field_Tester.cu
//---------------------------------------------------------------------------//
