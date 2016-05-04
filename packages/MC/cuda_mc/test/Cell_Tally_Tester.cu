//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Cell_Tally_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Cell_Tally_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Cell_Tally_Tester.hh"
#include "Test_XS.hh"
#include "../Cell_Tally.cuh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "xs/XS.hh"
#include "cuda_xs/XS_Device.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef profugus::XS                 XS_t;
typedef cuda_profugus::Mesh_Geometry Geom;
typedef cuda_profugus::Space_Vector  Space_Vector;

__global__ void test_tally_kernel( Cell_Tally<Geom> *tally,
                                   Geom             *geom,
                                   int               num_vals)
{
     using def::I; using def::J; using def::K;

     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         // Create particle
         Particle<Geom> p;
         p.set_group(0);
         p.set_matid(0);
         p.set_wt(1.0);
         p.live();

         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(tid,0,0,&rng_state);
         p.set_rng(&rng_state);
         
         // Initialize particle uniformly on [0, 1]
         double x_loc = curand_uniform_double(&rng_state);
         double y_loc = curand_uniform_double(&rng_state);
         double z_loc = curand_uniform_double(&rng_state);
         Space_Vector pos = {x_loc, y_loc, z_loc};

         // Direction doesn't matter
         Space_Vector dir = {1.0, 0.0, 0.0};

         geom->initialize(pos,dir,p.geo_state());

         // Tally with step length of 1.0
         tally->accumulate(1.0,p);

         // Move particle to new location
         pos[I] = curand_uniform_double(&rng_state);
         pos[J] = curand_uniform_double(&rng_state);
         pos[K] = curand_uniform_double(&rng_state);

         geom->initialize(pos,dir,p.geo_state());

         // Tally with step length of 0.5
         tally->accumulate(0.5,p);
     }
}

void Cell_Tally_Tester::test_tally()
{
    // Mesh edges
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};
    
    // Build geometry
    auto geom = std::make_shared<Geom>(x_edges,
        y_edges,z_edges);
    std::vector<int> matids(geom->num_cells(),0);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<Geom> sdp_geom(geom);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    auto xs = Test_XS::build_xs(1);

    int num_particles = 8192;
    std::vector<int> cells = {0, 1, 2, 4, 7};

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    auto sp_phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    sp_phys->set_geometry(sdp_geom);
    cuda::Shared_Device_Ptr<Physics<Geom> > sdp_phys(sp_phys);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Build cell tally
    auto sp_tally = std::make_shared<Cell_Tally<Geom> >(sdp_geom,sdp_phys);
    sp_tally->set_cells(cells);
    cuda::Shared_Device_Ptr<Cell_Tally<Geom> > tally(sp_tally);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Launch kernel
    int block_size = 256;
    int num_blocks = num_particles / block_size;
    if( num_blocks > 0 )
    {
        test_tally_kernel<<<num_blocks,block_size>>>(
            tally.get_device_ptr(), sdp_geom.get_device_ptr(), num_particles );
    }

    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaDeviceSynchronize();

    if( num_particles % block_size != 0 )
    {
        test_tally_kernel<<<1,num_particles%block_size>>>(
            tally.get_device_ptr(), sdp_geom.get_device_ptr(), num_particles );
    }

    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaDeviceSynchronize();

    tally.update_host();
    sp_tally->finalize(num_particles);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy tally result to host
    auto tally_result = sp_tally->results();
    EXPECT_EQ( tally_result.size(), cells.size() );

    // Each value should be 1.5 within statistical noise
    std::vector<double> expected(cells.size(),1.5);
    double tol = 10.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_VEC_SOFTEQ( expected, tally_result, tol );
    
}

//---------------------------------------------------------------------------//
//                 end of Cell_Tally_Tester.cc
//---------------------------------------------------------------------------//
