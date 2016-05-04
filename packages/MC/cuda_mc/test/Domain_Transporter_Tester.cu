//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Domain_Transporter_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Domain_Transporter_Tester.hh"
#include "../Domain_Transporter.cuh"
#include "../Uniform_Source.cuh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry      Geom;
typedef cuda_mc::Uniform_Source<Geom>     Uniform_Src;
typedef cuda_mc::Domain_Transporter<Geom> Transporter;
typedef profugus::XS                      XS_t;

__global__ void test_transport_kernel( Uniform_Src *source,
                                       Transporter *trans,
                                       int         *events,
                                       int          num_particles )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_particles )
     {
         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);

         // Get particle from source
         auto p = source->get_particle(tid,&rng_state);

         // Transport particle
         trans->transport(p);

         // Get final event
         events[tid] = p.event();
     }
}

namespace 
{

Teuchos::RCP<XS_t> build_5grp_xs()
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

Teuchos::RCP<XS_t> build_3grp_xs()
{
    const int ng = 3;
    Teuchos::RCP<XS_t> xs = Teuchos::rcp(new XS_t());
    xs->set(0, ng);

    // make group boundaries
    XS_t::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 1.0;
    nbnd[2] = 0.01;
    nbnd[3] = 0.0001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.1, 1.6, 2.9};
    double t2[ng] = {10.0, 11.3, 16.2};

    XS_t::OneDArray tot1(std::begin(t1), std::end(t1));
    XS_t::OneDArray tot2(std::begin(t2), std::end(t2));

    xs->add(0, XS_t::TOTAL, tot1);
    xs->add(1, XS_t::TOTAL, tot2);

    double s1[][3] = {{0.7, 0.0, 0.0},
                      {0.2, 0.3, 0.0},
                      {0.1, 0.7, 1.9}};

    double s2[][3] = {{2.7, 0.0, 0.0},
                      {2.2, 2.3, 0.0},
                      {2.1, 2.7, 3.9}};

    XS_t::TwoDArray sct1(ng, ng, 0.0);
    XS_t::TwoDArray sct2(ng, ng, 0.0);

    for (int g = 0; g < 3; ++g)
    {
        for (int gp = 0; gp < 3; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
            sct2(g, gp) = s2[g][gp];
        }
    }

    xs->add(0, 0, sct1);
    xs->add(1, 0, sct2);

    double c2[] = {0.4, 0.6, 0.0};
    double f2[] = {3.2, 4.2, 0.0};
    double n2[] = {2.4*3.2, 2.4*4.2, 0.0};

    XS_t::OneDArray chi2(std::begin(c2), std::end(c2));
    XS_t::OneDArray fis2(std::begin(f2), std::end(f2));
    XS_t::OneDArray nuf2(std::begin(n2), std::end(n2));

    xs->add(1, XS_t::CHI, chi2);
    xs->add(1, XS_t::SIG_F, fis2);
    xs->add(1, XS_t::NU_SIG_F, nuf2);

    xs->complete();

    return xs;
}

}
void Domain_Transporter_Tester::test_transport(int num_groups)
{
    REQUIRE(num_groups==3 || num_groups==5);

    int num_particles = 128;

    Teuchos::RCP<XS_t> xs;
    if (num_groups==3)
        xs = build_3grp_xs();
    else
        xs = build_5grp_xs();

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom = std::make_shared<Geom>(edges,edges,edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    pl->set("num_groups",xs->num_groups());
    pl->set("Np",num_particles);
    pl->set("implicit_capture",false);

    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);
    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    cuda::Shared_Device_Ptr<Physics<Geom> > sdp_phys(phys);

    // Build domain transporter
    auto transp = cuda::shared_device_ptr<Transporter>(sdp_geom,sdp_phys);

    // Build box shape for source
    std::vector<double> src_bounds = {edges.front(), edges.back(),
                                      edges.front(), edges.back(),
                                      edges.front(), edges.back()};
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            src_bounds[0], src_bounds[1],
            src_bounds[2], src_bounds[3],
            src_bounds[4], src_bounds[5]);

    // Build source
    auto sp_source = std::make_shared<Uniform_Src>(pl,sdp_geom);
    sp_source->build_source(src_shape);
    cuda::Shared_Device_Ptr<Uniform_Src> sdp_source(sp_source);

    // Allocate data on device
    thrust::device_vector<int> device_events(num_particles);

    test_transport_kernel<<<1,num_particles>>>( sdp_source.get_device_ptr(),
                                                transp.get_device_ptr(),
                                                device_events.data().get(),
                                                num_particles );

    cudaDeviceSynchronize();
    REQUIRE( cudaGetLastError() == cudaSuccess );

    thrust::host_vector<int> host_events = device_events;

    int num_absorptions = std::count_if(host_events.begin(),host_events.end(),
            [](int e){return e == profugus::events::ABSORPTION;});
    int num_escapes = std::count_if(host_events.begin(),host_events.end(),
            [](int e){return e == profugus::events::ESCAPE;});

    int expected_absorptions;
    int expected_escapes;
    // Heuristic tests...
    if (num_groups == 5)
    {
        expected_absorptions = 115;
        expected_escapes = 13;
    }
    else
    {
        expected_absorptions = 70;
        expected_escapes = 58;
    }
    EXPECT_EQ(expected_absorptions,num_absorptions);
    EXPECT_EQ(expected_escapes,num_escapes);
}

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter_Tester.cc

//---------------------------------------------------------------------------//
