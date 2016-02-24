//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Physics_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics_Tester.hh"
#include "../Physics.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry Geom;

__global__ void test_total_kernel( Physics<Geom> phys,
                                   double       *totals,
                                   int           num_vals)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         int g = tid % 5;
         int matid = tid % 2;

         // Create particle
         Particle<Geom> p;
         p.set_group(g);
         p.set_matid(matid);
         totals[tid] = phys.total(profugus::physics::TOTAL,p);
     }
}

__global__ void test_collide_kernel( Geom         *geom,
                                     Physics<Geom> phys,
                                     int           num_particles)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_particles )
     {
         int g = tid % 5;

         // Create particle
         Particle<Geom> p;
         p.live();
         p.set_group(g);
         p.set_wt(1.0);
         p.set_event(profugus::events::COLLISION);

         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);
         p.set_rng(rng_state);

         // Initialize geo state
         cuda::Space_Vector pos = {0.25, 0.75, 0.60};
         cuda::Space_Vector dir = {1.0, 0.0, 0.0};
         geom->initialize(pos,dir,p.geo_state());
         p.set_matid(geom->matid(p.geo_state()));

         // Collide
         phys.collide(p);

         printf("Particle %i has event %i, group %i, and weight %e\n",tid,p.event(),p.group(),p.wt());
     }
}

void Physics_Tester::test_total( const Vec_Dbl  &x_edges,
                                 const Vec_Dbl  &y_edges,
                                 const Vec_Dbl  &z_edges,
                                 const Vec_UInt &matids,
                                       SP_XS     xs,
                                       Vec_Dbl  &host_totals )
{

    // Build geometry
    auto geom = std::make_shared<cuda_profugus::Mesh_Geometry>(x_edges,
        y_edges,z_edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Physics<cuda_profugus::Mesh_Geometry> phys(pl,xs);
    phys.set_geometry(sdp_geom);

    // Allocate data on device
    int num_vals = host_totals.size();
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,double> device_totals(num_vals);

    test_total_kernel<<<1,num_vals>>>( phys, device_totals.data(),
                                       num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    device_totals.to_host(profugus::make_view(host_totals));
}

void Physics_Tester::test_collide( const Vec_Dbl  &x_edges,
                                   const Vec_Dbl  &y_edges,
                                   const Vec_Dbl  &z_edges,
                                   const Vec_UInt &matids,
                                         SP_XS     xs,
                                         int       num_particles )
{

    // Build geometry
    auto geom = std::make_shared<cuda_profugus::Mesh_Geometry>(x_edges,
        y_edges,z_edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Physics<cuda_profugus::Mesh_Geometry> phys(pl,xs);
    phys.set_geometry(sdp_geom);

    // Allocate data on device
    //typedef cuda::arch::Device Arch;
    //cuda::Device_Vector<Arch,double> device_totals(num_particles);

    test_collide_kernel<<<1,num_particles>>>( sdp_geom.get_device_ptr(),
                                              phys, num_particles );

    REQUIRE( cudaGetLastError() == cudaSuccess );
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Physics_Tester.cc
//---------------------------------------------------------------------------//
