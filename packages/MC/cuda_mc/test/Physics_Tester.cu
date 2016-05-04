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
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry Geom;
typedef cuda_profugus::Space_Vector  Space_Vector;
typedef profugus::XS                 XS_t;

__global__ void test_total_kernel( Physics<Geom> *phys,
                                   double        *totals,
                                   int            num_vals)
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
         totals[tid] = phys->total(profugus::physics::TOTAL,p);
     }
}

__global__ void test_collide_kernel( Geom          *geom,
                                     Physics<Geom> *phys,
                                     int           *events,
                                     int            num_particles)
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
         p.set_rng(&rng_state);

         // Initialize geo state
         Space_Vector pos = {0.25, 0.75, 0.60};
         Space_Vector dir = {1.0, 0.0, 0.0};
         geom->initialize(pos,dir,p.geo_state());
         p.set_matid(geom->matid(p.geo_state()));

         // Collide
         phys->collide(p);

         events[tid] = p.event();

         printf("Particle %i has event %i, group %i, and weight %e\n",
                tid,p.event(),p.group(),p.wt());
     }
}

namespace
{

Teuchos::RCP<XS_t> build_xs()
{
    Teuchos::RCP<XS_t> xs = Teuchos::rcp(new XS_t());
    xs->set(0, 5);

    std::vector<double> bnd(6, 0.0);
    bnd[0] = 100.0;
    bnd[1] = 10.0;
    bnd[2] = 1.0;
    bnd[3] = 0.1;
    bnd[4] = 0.01;
    bnd[5] = 0.001;
    xs->set_bounds(bnd);

    typename XS_t::OneDArray total(5);
    typename XS_t::TwoDArray scat(5, 5);

    double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
    double c[5] = {0.3770, 0.4421, 0.1809, 0.0, 0.0};
    double n[5] = {2.4*f[0], 2.4*f[1], 2.4*f[2], 2.4*f[3], 2.4*f[4]};
    typename XS_t::OneDArray fission(std::begin(f), std::end(f));
    typename XS_t::OneDArray chi(std::begin(c), std::end(c));
    typename XS_t::OneDArray nus(std::begin(n), std::end(n));
    xs->add(1, XS_t::SIG_F, fission);
    xs->add(1, XS_t::NU_SIG_F, nus);
    xs->add(1, XS_t::CHI, chi);

    // mat 0
    total[0] = 5.2 ;
    total[1] = 11.4;
    total[2] = 18.2;
    total[3] = 29.9;
    total[4] = 27.3;
    xs->add(0, XS_t::TOTAL, total);

    // mat 1
    total[0] = 5.2  + f[0];
    total[1] = 11.4 + f[1];
    total[2] = 18.2 + f[2];
    total[3] = 29.9 + f[3];
    total[4] = 27.3 + f[4];
    xs->add(1, XS_t::TOTAL, total);

    scat(0, 0) = 1.2;
    scat(1, 0) = 0.9;
    scat(1, 1) = 3.2;
    scat(2, 0) = 0.4;
    scat(2, 1) = 2.8;
    scat(2, 2) = 6.9;
    scat(2, 3) = 1.5;
    scat(3, 0) = 0.1;
    scat(3, 1) = 2.1;
    scat(3, 2) = 5.5;
    scat(3, 3) = 9.7;
    scat(3, 4) = 2.1;
    scat(4, 1) = 0.2;
    scat(4, 2) = 1.3;
    scat(4, 3) = 6.6;
    scat(4, 4) = 9.9;
    xs->add(0, 0, scat);
    xs->add(1, 0, scat);

    xs->complete();

    return xs;
}

}

void Physics_Tester::test_total()
{
    auto xs = build_xs();

    int num_vals = 16;

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom = std::make_shared<cuda_profugus::Mesh_Geometry>(edges,
        edges,edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);
    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    auto sdp_phys = cuda::Shared_Device_Ptr<Physics<Geom> >(phys);

    // Allocate data on device
    thrust::device_vector<double> device_totals(num_vals);

    test_total_kernel<<<1,num_vals>>>( sdp_phys.get_device_ptr(),
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
    auto xs = build_xs();

    int num_particles = 16;

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom = std::make_shared<cuda_profugus::Mesh_Geometry>(edges,
        edges,edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);

    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    auto sdp_phys = cuda::Shared_Device_Ptr<Physics<Geom> >(phys);

    thrust::device_vector<int> device_events(num_particles);

    test_collide_kernel<<<1,num_particles>>>( sdp_geom.get_device_ptr(),
                                              sdp_phys.get_device_ptr(),
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
