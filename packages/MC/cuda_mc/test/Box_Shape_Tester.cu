//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Box_Shape_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Box_Shape_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Box_Shape_Tester.hh"
#include "../Box_Shape.cuh"

#include "utils/View_Field.hh"
#include "cuda_utils/Device_Vector.hh"

namespace cuda_mc
{

__global__ void compute_inside_kernel( Box_Shape                 box,
                                       const cuda::Space_Vector *pts,
                                       int                      *inside,
                                       int                       num_vals )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         inside[tid] = box.is_point_inside(pts[tid]);
     }
}

__global__ void compute_sample_kernel( Box_Shape           box,
                                       cuda::Space_Vector *pts,
                                       int                 num_vals )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);
         pts[tid] = box.sample(rng_state);
     }
}

void Box_Shape_Tester::test_inside( const Vec_Dbl       &box_bounds,
                                    const Vec_Space_Vec &pts_host,
                                    Vec_Int             &inside_host)
{
    // Copy values to device
    int num_vals = pts_host.size();
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cuda::Space_Vector> pts_device(num_vals);
    pts_device.assign(profugus::make_view(pts_host));

    // Build box to be copied to device
    REQUIRE( box_bounds.size() == 6 );
    cuda_mc::Box_Shape box = {box_bounds[0], box_bounds[1],
                              box_bounds[2], box_bounds[3],
                              box_bounds[4], box_bounds[5]};

    // Storage for result of kernel
    cuda::Device_Vector<Arch,int> inside_device(num_vals);

    // Launch kernel
    compute_inside_kernel<<<1,num_vals>>>( box, pts_device.data(),
                                           inside_device.data(), num_vals );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    inside_device.to_host(profugus::make_view(inside_host));
}

void Box_Shape_Tester::test_sample( const Vec_Dbl &box_bounds,
                                    Vec_Space_Vec &pts_host)
{
    // Copy values to device
    int num_vals = pts_host.size();
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cuda::Space_Vector> pts_device(num_vals);

    // Build box to be copied to device
    REQUIRE( box_bounds.size() == 6 );
    cuda_mc::Box_Shape box = {box_bounds[0], box_bounds[1],
                              box_bounds[2], box_bounds[3],
                              box_bounds[4], box_bounds[5]};

    // Launch kernel
    compute_sample_kernel<<<1,num_vals>>>( box, pts_device.data(), num_vals );
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    pts_device.to_host(profugus::make_view(pts_host));
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Box_Shape_Tester.cc
//---------------------------------------------------------------------------//
