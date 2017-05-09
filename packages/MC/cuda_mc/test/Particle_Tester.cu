//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Particle_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Particle_Tester.hh"
#include "../Particle_Vector.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "gtest/Gtest_Functions.hh"

using namespace cuda_mc;

__global__ void compute_randoms_kernel(
        Particle_Vector<cuda_profugus::Mesh_Geometry> particles,
        double *randoms, int num_vals )

{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         // Generate random
         randoms[tid] = particles.ran(tid);
     }
}

__global__ void compute_groups_kernel( int *groups_in,
                                       int *groups_out,
                                       int num_vals )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         // Create particle and assign group
         Particle<cuda_profugus::Mesh_Geometry> p;
         p.set_group(groups_in[tid]);

         // Generate random
         groups_out[tid] = p.group();
     }
}

__global__ void compute_matids_kernel( int *matids_in,
        int *matids_out,
        int num_vals )
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_vals )
    {
        // Create particle and assign matid
        Particle<cuda_profugus::Mesh_Geometry> p;
        p.set_matid(matids_in[tid]);

        // Generate random
        matids_out[tid] = p.matid();
    }
}

void Particle_Tester::test_randoms()
{
    for (int num_vals : {16, 64, 256, 1024})
    {
        thrust::device_vector<double> device_rands(num_vals);

        Particle_Vector_DMM<cuda_profugus::Mesh_Geometry> particles(1234);

        compute_randoms_kernel<<<1,num_vals>>>(particles.device_instance(),
                                               device_rands.data().get(),
                                               num_vals );

        REQUIRE( cudaGetLastError() == cudaSuccess );

        thrust::host_vector<double> host_rands = device_rands;

        // Test statistics on values
        double mean = 0.0;
        for (auto val : host_rands)
        {
            EXPECT_GT(val,0.0);
            EXPECT_LT(val,1.0);
            mean += val;
        }
        double N = static_cast<double>(num_vals);
        mean /= N;

        EXPECT_SOFTEQ(0.5,mean,1.0/std::sqrt(N));
    }
}

void Particle_Tester::test_groups()
{
    int num_vals = 16;

    thrust::host_vector<int> host_groups_in(num_vals);
    for (int i = 0; i < num_vals; ++i)
        host_groups_in[i] = i % 4;

    thrust::device_vector<int> device_groups_in = host_groups_in;
    thrust::device_vector<int> device_groups_out(num_vals);

    compute_groups_kernel<<<1,num_vals>>>( device_groups_in.data().get(), 
                                           device_groups_out.data().get(),
                                           num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    thrust::host_vector<int> host_groups_out = device_groups_out;
    for (int i = 0; i < num_vals; ++i)
        EXPECT_EQ(i%4,host_groups_out[i]);
}

void Particle_Tester::test_matids()
{
    int num_vals = 16;

    thrust::host_vector<int> host_matids_in(num_vals);
    for (int i = 0; i < num_vals; ++i)
        host_matids_in[i] = i % 4;

    thrust::device_vector<int> device_matids_in = host_matids_in;
    thrust::device_vector<int> device_matids_out(num_vals);

    compute_matids_kernel<<<1,num_vals>>>( device_matids_in.data().get(), 
                                           device_matids_out.data().get(),
                                           num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    thrust::host_vector<int> host_matids_out = device_matids_out;
    for (int i = 0; i < num_vals; ++i)
        EXPECT_EQ(i%4,host_matids_out[i]);
}

//---------------------------------------------------------------------------//
//                 end of Particle_Tester.cc
//---------------------------------------------------------------------------//
