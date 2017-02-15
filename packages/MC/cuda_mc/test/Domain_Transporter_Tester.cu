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
#include "Test_XS.hh"
#include "../Domain_Transporter.cuh"
#include "../Uniform_Source.cuh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry      Geom;
typedef cuda_profugus::Mesh_Geometry_DMM  Geom_DMM;
typedef cuda_mc::Uniform_Source<Geom>     Uniform_Src;
typedef cuda_mc::Uniform_Source_DMM<Geom> Uniform_Src_DMM;
typedef cuda_mc::Domain_Transporter<Geom> Transporter;
typedef cuda_mc::Domain_Transporter_DMM<Geom> Transporter_DMM;

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

void Domain_Transporter_Tester::test_transport(int num_groups)
{
    REQUIRE(num_groups==3 || num_groups==5);

    def::size_type num_particles = 128;

    auto xs = Test_XS::build_xs(num_groups);

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom_dmm = std::make_shared<Geom_DMM>(edges,edges,edges);
    geom_dmm->set_matids(matids);
    auto sdp_geom = cuda::shared_device_ptr<Geom>(geom_dmm->device_instance());

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
    auto transp_dmm = std::make_shared<Transporter_DMM>(pl,sdp_geom,sdp_phys);
    auto transp = cuda::shared_device_ptr<Transporter>(transp_dmm->device_instance());

    // Build box shape for source
    std::vector<double> src_bounds = {edges.front(), edges.back(),
                                      edges.front(), edges.back(),
                                      edges.front(), edges.back()};
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            src_bounds[0], src_bounds[1],
            src_bounds[2], src_bounds[3],
            src_bounds[4], src_bounds[5]);

    // Build source
    auto sp_source = std::make_shared<Uniform_Src_DMM>(pl,sdp_geom);
    sp_source->build_source(src_shape);
    auto sdp_source = cuda::shared_device_ptr<Uniform_Src>(
            sp_source->device_instance());

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
