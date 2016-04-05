//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Keff_Tally_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Keff_Tally_Tester.hh"
#include "../Keff_Tally.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry Geom;

__global__ void test_tally_kernel( Keff_Tally<Geom> *tally,
                                   Geom             *geom,
                                   int               num_groups,
                                   int               num_vals)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         // Create particle
         Particle<Geom> p;
         p.set_group(tid % num_groups);
         p.set_matid(0);
         p.set_wt(1.0);

         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);
         p.set_rng(&rng_state);
         
         // Initialize particle uniformly on [0, 1]
         double x_loc = curand_uniform_double(&rng_state);
         double y_loc = curand_uniform_double(&rng_state);
         double z_loc = curand_uniform_double(&rng_state);
         cuda::Space_Vector pos = {x_loc, y_loc, z_loc};

         // Direction doesn't matter
         cuda::Space_Vector dir = {1.0, 0.0, 0.0};

         geom->initialize(pos,dir,p.geo_state());

         // Tally with step length of 1.0
         tally->accumulate(1.0,p);

         // Move particle to new location
         pos.x = curand_uniform_double(&rng_state);
         pos.y = curand_uniform_double(&rng_state);
         pos.z = curand_uniform_double(&rng_state);

         geom->initialize(pos,dir,p.geo_state());

         // Tally with step length of 0.5
         tally->accumulate(0.5,p);
     }
}

void Keff_Tally_Tester::test_tally( const Vec_Dbl  &x_edges,
                                    const Vec_Dbl  &y_edges,
                                    const Vec_Dbl  &z_edges,
                                          RCP_XS     xs,
                                          double   &keff,
                                          int       num_particles )
{
    // Build geometry
    auto geom = std::make_shared<Geom>(x_edges,
        y_edges,z_edges);
    std::vector<int> matids(geom->num_cells(),0);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<Geom> sdp_geom(geom);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    auto sp_phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    sp_phys->set_geometry(sdp_geom);
    cuda::Shared_Device_Ptr<Physics<Geom> > sdp_phys(sp_phys);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Build cell tally
    std::cout << "Building Keff_Tally" << std::endl;
    auto sp_tally = std::make_shared<Keff_Tally<Geom> >(1.0,sdp_phys);
    cuda::Shared_Device_Ptr<Keff_Tally<Geom> > tally(sp_tally);

    sp_tally->begin_active_cycles();
    sp_tally->begin_cycle(tally.get_device_ptr());

    // Launch kernel
    std::cout << "Launching kernel" << std::endl;
    int block_size = 256;
    int num_blocks = num_particles / block_size;
    if( num_blocks > 0 )
    {
        test_tally_kernel<<<num_blocks,block_size>>>( tally.get_device_ptr(),
                                                      sdp_geom.get_device_ptr(),
                                                      sp_phys->num_groups(),
                                                      num_particles );
    }

    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaDeviceSynchronize();

    if( num_particles % block_size != 0 )
    {
        test_tally_kernel<<<1,num_particles%block_size>>>( tally.get_device_ptr(),
                                                           sdp_geom.get_device_ptr(),
                                                           sp_phys->num_groups(),
                                                           num_particles );
    }

    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaDeviceSynchronize();

    sp_tally->end_cycle(num_particles,tally.get_device_ptr());

    // Copy tally result to host
    keff = sp_tally->mean();
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Keff_Tally_Tester.cc
//---------------------------------------------------------------------------//
