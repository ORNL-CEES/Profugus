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
#include "../Physics.hh"
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

void Physics_Tester::test_total( const Vec_Dbl  &x_edges,
                                 const Vec_Dbl  &y_edges,
                                 const Vec_Dbl  &z_edges,
                                 const Vec_UInt &matids,
                                       SP_XS     xs,
                                       Vec_Dbl  &host_totals )
{

    // Build geometry
    cuda_profugus::Mesh_Geometry geom(x_edges,y_edges,z_edges);
    geom.set_matids(matids);
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

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Physics_Tester.cc
//---------------------------------------------------------------------------//
