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
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry      Geom;
typedef cuda_mc::Uniform_Source<Geom>     Uniform_Src;
typedef cuda_mc::Domain_Transporter<Geom> Transporter;

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
         auto p = source->get_particle(tid,rng_state);

         // Transport particle
         trans->transport(p);

         // Get final event
         events[tid] = p.event();
     }
}

void Domain_Transporter_Tester::test_transport( const Vec_Dbl  &x_edges,
                                                const Vec_Dbl  &y_edges,
                                                const Vec_Dbl  &z_edges,
                                                const Vec_UInt &matids,
                                                      SP_XS     xs,
                                                      Vec_Int  &events )
{
    int num_particles = events.size();

    // Build geometry
    auto geom = std::make_shared<Geom>(x_edges,y_edges,z_edges);
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
    auto sp_trans = std::make_shared<Transporter>();
    sp_trans->set(sdp_geom,sdp_phys);
    cuda::Shared_Device_Ptr<Transporter> sdp_trans(sp_trans);

    // Build box shape for source
    Vec_Dbl src_bounds = {x_edges.front(), x_edges.back(),
                          y_edges.front(), y_edges.back(),
                          z_edges.front(), z_edges.back()};
    REQUIRE( src_bounds.size() == 6 );
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            src_bounds[0], src_bounds[1],
            src_bounds[2], src_bounds[3],
            src_bounds[4], src_bounds[5]);

    // Build source
    auto sp_source = std::make_shared<Uniform_Src>(pl,sdp_geom);
    sp_source->build_source(src_shape);
    cuda::Shared_Device_Ptr<Uniform_Src> sdp_source(sp_source);

    // Allocate data on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int> device_events(num_particles);

    test_transport_kernel<<<1,num_particles>>>( sdp_source.get_device_ptr(),
                                                sdp_trans.get_device_ptr(),
                                                device_events.data(),
                                                num_particles );

    cudaDeviceSynchronize();
    REQUIRE( cudaGetLastError() == cudaSuccess );

    device_events.to_host(profugus::make_view(events));
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter_Tester.cc

//---------------------------------------------------------------------------//
