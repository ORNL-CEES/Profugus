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
#include "../Particle.cuh"
#include "../../cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

 __global__ void compute_randoms_kernel( double *randoms, int num_vals )
 {
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);

         // Create particle and assign RNG
         Particle<cuda_profugus::Mesh_Geometry> p;
         p.set_rng(rng_state);

         // Generate random
         randoms[tid] = p.ran();
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

void Particle_Tester::test_randoms( Vec_Dbl &host_rands )
{
    int num_vals = host_rands.size();

    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,double> device_rands(num_vals);

    compute_randoms_kernel<<<1,num_vals>>>( device_rands.data(), num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    device_rands.to_host(profugus::make_view(host_rands));
}

void Particle_Tester::test_groups( const Vec_Int &host_groups_in,
                                         Vec_Int &host_groups_out)
{
    int num_vals = host_groups_in.size();
    REQUIRE( host_groups_out.size() == num_vals );

    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int> device_groups_in(num_vals);
    cuda::Device_Vector<Arch,int> device_groups_out(num_vals);

    device_groups_in.assign(profugus::make_view(host_groups_in));

    compute_groups_kernel<<<1,num_vals>>>( device_groups_in.data(), 
                                           device_groups_out.data(),
                                           num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    device_groups_out.to_host(profugus::make_view(host_groups_out));
}

void Particle_Tester::test_matids( const Vec_Int &host_matids_in,
                                         Vec_Int &host_matids_out)
{
    int num_vals = host_matids_in.size();
    REQUIRE( host_matids_out.size() == num_vals );

    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int> device_matids_in(num_vals);
    cuda::Device_Vector<Arch,int> device_matids_out(num_vals);

    device_matids_in.assign(profugus::make_view(host_matids_in));

    compute_matids_kernel<<<1,num_vals>>>( device_matids_in.data(), 
                                           device_matids_out.data(),
                                           num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    device_matids_out.to_host(profugus::make_view(host_matids_out));
}
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Particle_Tester.cc
//---------------------------------------------------------------------------//
