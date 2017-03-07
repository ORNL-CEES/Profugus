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

typedef profugus::XS                        XS_t;
typedef cuda_profugus::Mesh_Geometry        Geom;
typedef cuda_profugus::Mesh_Geometry_DMM    Geom_DMM;
typedef cuda_utils::Space_Vector            Space_Vector;
typedef Particle_Vector<Geom>               Particle_Vector_t;
typedef Particle_Vector_DMM<Geom>           Particle_Vector_DMM_t;

__global__ void test_tally_kernel( Cell_Tally<Geom> *tally,
                                   Geom             *geom,
                                   Particle_Vector_t particles,
                                   int               num_vals,
                                   int               seed)
{
     using def::I; using def::J; using def::K;

     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         // Create particle
         particles.set_group(tid,0);
         particles.set_matid(tid,0);
         particles.set_wt(tid,1.0);
         particles.live(tid);

         // Create and initialize RNG state
         curandState_t rng_state;
         curand_init(seed,tid,0,&rng_state);
         particles.set_rng(tid,&rng_state);
         
         // Initialize particle uniformly on [0, 1]
         double x_loc = curand_uniform_double(&rng_state);
         double y_loc = curand_uniform_double(&rng_state);
         double z_loc = curand_uniform_double(&rng_state);
         Space_Vector pos = {x_loc, y_loc, z_loc};

         // Direction doesn't matter
         Space_Vector dir = {1.0, 0.0, 0.0};

         geom->initialize(pos,dir,particles.geo_state(tid));

         // Tally with step length of 1.0
         tally->accumulate(1.0,tid,particles);

         // Move particle to new location
         pos[I] = curand_uniform_double(&rng_state);
         pos[J] = curand_uniform_double(&rng_state);
         pos[K] = curand_uniform_double(&rng_state);

         geom->initialize(pos,dir,particles.geo_state(tid));

         // Tally with step length of 0.5
         tally->accumulate(0.5,tid,particles);
     }
}

void Cell_Tally_Tester::test_tally(int num_batches)
{
    // Mesh edges
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};
    
    // Build geometry
    auto geom_dmm = std::make_shared<Geom_DMM>(x_edges,
        y_edges,z_edges);
    std::vector<int> matids(geom_dmm->num_cells(),0);
    geom_dmm->set_matids(matids);
    auto sdp_geom = cuda::shared_device_ptr<Geom>(geom_dmm->device_instance());

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
    auto sp_tally_dmm = std::make_shared<Cell_Tally_DMM<Geom> >(sdp_geom,sdp_phys);
    sp_tally_dmm->set_cells(cells,geom_dmm->volumes());
    auto tally = cuda::shared_device_ptr<Cell_Tally<Geom>>(
        sp_tally_dmm->device_instance());

    // build particles
    Particle_Vector_DMM_t particles;
    particles.initialize(num_particles);

    REQUIRE( cudaGetLastError() == cudaSuccess );

    for (int batch = 0; batch < num_batches; ++batch)
    {
        // Launch kernel
        int block_size = 256;
        int num_blocks = num_particles / block_size;
        REQUIRE(num_particles % block_size == 0);
        if( num_blocks > 0 )
        {
            test_tally_kernel<<<num_blocks,block_size>>>(
                tally.get_device_ptr(), sdp_geom.get_device_ptr(),
                particles.device_instance(),
                num_particles, batch);
        }

        REQUIRE( cudaGetLastError() == cudaSuccess );
        cudaDeviceSynchronize();

        sp_tally_dmm->end_batch(num_particles);
    }

    sp_tally_dmm->finalize(num_particles*num_batches);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Get tally result
    const auto &tally_result = sp_tally_dmm->results();
    EXPECT_EQ( tally_result.size(), cells.size() );

    const auto &tally_std_dev = sp_tally_dmm->std_dev();
    EXPECT_EQ( tally_std_dev.size(), cells.size() );

    // Each value should be 1.5 within statistical noise
    std::vector<double> expected(cells.size(),1.5);

    // If only one batch then the computed variance is zero,
    //  evaluate heuristically
    if (num_batches == 1)
    {
        double tol = 10.0 / std::sqrt( static_cast<double>(num_particles*num_batches) );
        EXPECT_VEC_SOFTEQ( expected, tally_result, tol );

        std::vector<double> z(cells.size(),0.0);
        EXPECT_VEC_SOFT_EQ(z, tally_std_dev);
    }
    else
    {
        // Make sure all values are within 3 sigma, half within 1 sigma
        int num_converged = 0;
        for (int cell = 0; cell < cells.size(); ++cell)
        {
            double sig = tally_std_dev[cell];
            EXPECT_SOFTEQ(expected[cell],tally_result[cell],3*sig);
            if (profugus::soft_equiv(expected[cell],tally_result[cell],sig))
                num_converged++;
        }
        std::cout << num_converged << " out of " << cells.size()
            << " cells were within 1 sigma of true mean" << std::endl;
        EXPECT_GE(num_converged,cells.size()/2);
    }
}

//---------------------------------------------------------------------------//
//                 end of Cell_Tally_Tester.cc
//---------------------------------------------------------------------------//
