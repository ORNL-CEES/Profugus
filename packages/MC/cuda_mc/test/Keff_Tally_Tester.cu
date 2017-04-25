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
#include "Test_XS.hh"
#include "../Keff_Tally.cuh"
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

__global__ void test_tally_kernel( Keff_Tally<Geom> *tally,
                                   Geom             *geom,
                                   Particle_Vector_t particles,
                                   int               num_groups,
                                   int               num_vals)
{
    using def::I; using def::J; using def::K;

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_vals )
    {
        // Create particle
        particles.set_group(tid,tid % num_groups);
        particles.set_matid(tid,0);
        particles.set_wt(tid,1.0);

        // Create and initialize RNG state
        curandState_t rng_state;
        curand_init(tid,0,0,&rng_state);
        particles.set_rng(tid,&rng_state);

        // Initialize particle uniformly on [0, 1]
        double x_loc = curand_uniform_double(&rng_state);
        double y_loc = curand_uniform_double(&rng_state);
        double z_loc = curand_uniform_double(&rng_state);
        Space_Vector pos = {x_loc, y_loc, z_loc};

        // Direction doesn't matter
        Space_Vector dir = {1.0, 0.0, 0.0};

        geom->initialize(pos,dir,particles.geo_states(),tid);

        // Tally with step length of 1.0
        tally->accumulate(1.0,tid,particles);

        // Move particle to new location
        pos[I] = curand_uniform_double(&rng_state);
        pos[J] = curand_uniform_double(&rng_state);
        pos[K] = curand_uniform_double(&rng_state);

        geom->initialize(pos,dir,particles.geo_states(),tid);

        // Tally with step length of 0.5
        tally->accumulate(0.5,tid,particles);
    }
}

void Keff_Tally_Tester::test_tally(int num_groups)
{
    int num_particles = 64;

    REQUIRE( num_groups==1 || num_groups==2 );

    auto xs = Test_XS::build_xs(num_groups);
    REQUIRE(xs != Teuchos::null);

    // Build geometry
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};
    auto geom_dmm = std::make_shared<Geom_DMM>(x_edges,
        y_edges,z_edges);
    std::vector<int> matids(geom_dmm->num_cells(),0);
    geom_dmm->set_matids(matids);
    auto sdp_geom = cuda::shared_device_ptr<Geom>(geom_dmm->device_instance());

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
    auto sp_tally = std::make_shared<Keff_Tally_DMM<Geom> >(1.0,sdp_phys);

    sp_tally->begin_active_cycles();
    sp_tally->begin_cycle(num_particles);
    auto tally = cuda::shared_device_ptr<Keff_Tally<Geom>>(
        sp_tally->device_instance());

    // Build particles
    Particle_Vector_DMM_t particles;
    particles.initialize(num_particles);

    // Launch kernel
    std::cout << "Launching kernel" << std::endl;
    int block_size = 256;
    int num_blocks = num_particles / block_size;
    if( num_blocks > 0 )
    {
        test_tally_kernel<<<num_blocks,block_size>>>( tally.get_device_ptr(),
                                                      sdp_geom.get_device_ptr(),
                                                      particles.device_instance(),
                                                      sp_phys->num_groups(),
                                                      num_particles );
    }

    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaDeviceSynchronize();

    if( num_particles % block_size != 0 )
    {
        test_tally_kernel<<<1,num_particles%block_size>>>( tally.get_device_ptr(),
                                                           sdp_geom.get_device_ptr(),
                                                           particles.device_instance(),
                                                           sp_phys->num_groups(),
                                                           num_particles );
    }

    REQUIRE( cudaGetLastError() == cudaSuccess );
    cudaDeviceSynchronize();

    sp_tally->end_cycle(num_particles);

    // Copy tally result to host
    double keff = sp_tally->mean();

    double expected_keff;
    if (num_groups == 1)
    {
        // nu_sigma_f is 0.4 and each thread adds contribution of 1.5
        expected_keff = 1.5 * 0.4;
    }
    else
    {
        // Value should be 0.6
        // Half threads in group 0 add 1.5 * 0.4
        // Half threads in group 1 add 1.5 * 0.05
        expected_keff = 0.5 * (1.5 * 0.4 + 1.5 * 0.05);
    }
    double tol = 1.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_SOFTEQ( expected_keff, keff, tol );
        
}

//---------------------------------------------------------------------------//
//                 end of Keff_Tally_Tester.cc
//---------------------------------------------------------------------------//
