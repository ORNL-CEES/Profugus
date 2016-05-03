//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Box_Shape_Tester.cu
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Box_Shape_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>

#include "Box_Shape_Tester.hh"
#include "../Box_Shape.cuh"

#include "utils/View_Field.hh"
#include "gtest/Gtest_Functions.hh"

using namespace cuda_mc;

typedef cuda_utils::Space_Vector Space_Vector;

__global__ void compute_inside_kernel( Box_Shape           box,
                                       const Space_Vector *pts,
                                       int                *inside,
                                       int                 num_vals )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         inside[tid] = box.is_point_inside(pts[tid]);
     }
}

__global__ void compute_sample_kernel( Box_Shape     box,
                                       Space_Vector *pts,
                                       int           num_vals )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);
         pts[tid] = box.sample(&rng_state);
     }
}

void Box_Shape_Tester::test_inside()
{
    int num_vals = 4;
    std::vector<cuda_utils::Space_Vector> host_pts =
        {{1.2,   0.0, 3.5},
         {0.5,   0.1, 4.9},
         {0.75, -2.0, 4.0},
         {0.1,  -0.4, 5.01}};
    
    // Copy values to device
    thrust::device_vector<Space_Vector> device_pts(host_pts);

    // Build box to be copied to device
    cuda_mc::Box_Shape box = {0.0, 1.0, -1.0, 1.0, 3.0, 5.0};

    // Storage for result of kernel
    thrust::device_vector<int> device_inside(num_vals);

    // Launch kernel
    compute_inside_kernel<<<1,num_vals>>>( box,
                                           device_pts.data().get(),
                                           device_inside.data().get(),
                                           num_vals );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    std::vector<int> host_inside(num_vals);
    thrust::copy(device_inside.begin(),device_inside.end(),host_inside.begin());

    std::vector<int> expected_inside = {0, 1, 0, 0};
    EXPECT_VEC_EQ(expected_inside,host_inside);
}

void Box_Shape_Tester::test_sample()
{
    // Copy values to device
    int num_vals = 100;
    thrust::device_vector<Space_Vector> device_pts(num_vals);

    // Build box to be copied to device
    cuda_mc::Box_Shape box = {0.0, 1.0, -1.0, 1.0, 3.0, 5.0};

    // Launch kernel
    compute_sample_kernel<<<1,num_vals>>>(box, device_pts.data().get(),
                                          num_vals );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    std::vector<cuda_utils::Space_Vector> host_pts(num_vals);
    thrust::copy(device_pts.begin(),device_pts.end(),host_pts.begin());

    // Make sure all points are in bounding box
    for (const auto &pt : host_pts)
    {
        EXPECT_GE(pt.x,0.0);
        EXPECT_LE(pt.x,1.0);
        EXPECT_GE(pt.y,-1.0);
        EXPECT_LE(pt.y,1.0);
        EXPECT_GE(pt.z,3.0);
        EXPECT_LE(pt.z,5.0);
    }
}

//---------------------------------------------------------------------------//
//                 end of Box_Shape_Tester.cu
//---------------------------------------------------------------------------//
