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
#include "Test_XS.hh"
#include "../Physics.cuh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry     Geom;
typedef cuda_profugus::Mesh_Geometry_DMM Geom_DMM;
typedef cuda_utils::Space_Vector         Space_Vector;
typedef Particle_Vector<Geom>            Particle_Vector_t;
typedef Particle_Vector_DMM<Geom>        Particle_Vector_DMM_t;

__global__ void test_total_kernel( Physics<Geom>    *phys,
                                   Particle_Vector_t particles,
                                   double           *totals,
                                   int               num_vals)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         int g = tid % 5;
         int matid = tid % 2;

         // Create particle
         particles.set_group(tid,g);
         particles.set_matid(tid,matid);
         totals[tid] = phys->total(profugus::physics::TOTAL,tid,particles);
     }
}

__global__ void test_collide_kernel( Geom               *geom,
                                     Physics<Geom>      *phys,
                                     Particle_Vector_t   particles,
                                     int                *events,
                                     int                 num_particles)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_particles )
     {
         int g = tid % 5;

         // Create particle
         particles.live(tid);
         particles.set_group(tid,g);
         particles.set_wt(tid,1.0);
         particles.set_event(tid,profugus::events::COLLISION);

         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);
         particles.set_rng(tid,&rng_state);

         // Initialize geo state
         Space_Vector pos = {0.25, 0.75, 0.60};
         Space_Vector dir = {1.0, 0.0, 0.0};
         geom->initialize(pos,dir,particles.geo_state(tid));
         particles.set_matid(tid,geom->matid(particles.geo_state(tid)));

         // Collide
         phys->collide(tid,particles);

         events[tid] = particles.event(tid);

         printf("Particle %i has event %i, group %i, and weight %e\n",
            tid,particles.event(tid),particles.group(tid),particles.wt(tid));
     }
}

void Physics_Tester::test_total()
{
    auto xs = Test_XS::build_xs(5);

    int num_vals = 16;

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom_dmm = std::make_shared<Geom_DMM>(edges,edges,edges);
    geom_dmm->set_matids(matids);
    auto sdp_geom = cuda::shared_device_ptr<Geom>(geom_dmm->device_instance());

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);
    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    auto sdp_phys = cuda::Shared_Device_Ptr<Physics<Geom> >(phys);

    // Build particles
    Particle_Vector_DMM_t particles;
    particles.initialize(num_vals);

    // Allocate data on device
    thrust::device_vector<double> device_totals(num_vals);

    test_total_kernel<<<1,num_vals>>>( sdp_phys.get_device_ptr(),
                                       particles.device_instance(),
                                       device_totals.data().get(),
                                       num_vals );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    thrust::host_vector<double> host_totals = device_totals;

    for (int i = 0; i < num_vals; ++i)
    {
        int g = i % 5;
        int matid = i %2;
        const auto &expected = xs->vector(matid,profugus::XS::TOTAL);
        EXPECT_SOFT_EQ(expected[g], host_totals[i]);
    }
}

void Physics_Tester::test_collide()
{
    auto xs = Test_XS::build_xs(5);

    int num_particles = 16;

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom_dmm = std::make_shared<Geom_DMM>(edges,edges,edges);
    geom_dmm->set_matids(matids);
    auto sdp_geom = cuda::shared_device_ptr<Geom>(geom_dmm->device_instance());

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);

    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    auto sdp_phys = cuda::Shared_Device_Ptr<Physics<Geom> >(phys);

    // Build particles
    Particle_Vector_DMM_t particles;
    particles.initialize(num_particles);

    thrust::device_vector<int> device_events(num_particles);

    test_collide_kernel<<<1,num_particles>>>( sdp_geom.get_device_ptr(),
                                              sdp_phys.get_device_ptr(),
                                              particles.device_instance(),
                                              device_events.data().get(),
                                              num_particles );

    EXPECT_EQ(cudaGetLastError(), cudaSuccess);

    thrust::host_vector<int> host_events = device_events;
    for (auto event : host_events)
        EXPECT_EQ(profugus::events::IMPLICIT_CAPTURE,event);

}

//---------------------------------------------------------------------------//
//                 end of Physics_Tester.cc
//---------------------------------------------------------------------------//
