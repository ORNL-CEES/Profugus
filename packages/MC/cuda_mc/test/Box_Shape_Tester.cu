//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Box_Shape_Tester.cc
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

namespace cuda_mc
{

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

void Box_Shape_Tester::test_inside( const Vec_Dbl       &box_bounds,
                                    const Vec_Space_Vec &host_pts,
                                    Vec_Int             &host_inside)
{
    // Copy values to device
    int num_vals = host_pts.size();
    thrust::device_vector<Space_Vector> device_pts(host_pts);

    // Build box to be copied to device
    REQUIRE( box_bounds.size() == 6 );
    cuda_mc::Box_Shape box = {box_bounds[0], box_bounds[1],
                              box_bounds[2], box_bounds[3],
                              box_bounds[4], box_bounds[5]};

    // Storage for result of kernel
    thrust::device_vector<int> device_inside(num_vals);

    // Launch kernel
    compute_inside_kernel<<<1,num_vals>>>( box,
                                           device_pts.data().get(),
                                           device_inside.data().get(),
                                           num_vals );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    thrust::copy(device_inside.begin(),device_inside.end(),host_inside.begin());
}

void Box_Shape_Tester::test_sample( const Vec_Dbl &box_bounds,
                                    Vec_Space_Vec &host_pts)
{
    // Copy values to device
    int num_vals = host_pts.size();
    thrust::device_vector<Space_Vector> device_pts(num_vals);

    // Build box to be copied to device
    REQUIRE( box_bounds.size() == 6 );
    cuda_mc::Box_Shape box = {box_bounds[0], box_bounds[1],
                              box_bounds[2], box_bounds[3],
                              box_bounds[4], box_bounds[5]};

    // Launch kernel
    compute_sample_kernel<<<1,num_vals>>>(box, device_pts.data().get(),
                                          num_vals );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    thrust::copy(device_pts.begin(),device_pts.end(),host_pts.begin());
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Box_Shape_Tester.cc
//---------------------------------------------------------------------------//
