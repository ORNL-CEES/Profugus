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
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry Geom;
typedef cuda_profugus::Space_Vector  Space_Vector;
typedef profugus::XS                 XS_t;

__global__ void test_tally_kernel( Keff_Tally<Geom> *tally,
                                   Geom             *geom,
                                   int               num_groups,
                                   int               num_vals)
{
    using def::I; using def::J; using def::K;

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

namespace
{

Teuchos::RCP<XS_t> build_1grp_xs()
{
    constexpr int ng = 1;
    Teuchos::RCP<XS_t> xs = Teuchos::rcp(new XS_t());
    xs->set(0, ng);

    // make group boundaries
    XS_t::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 0.00001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.0};

    XS_t::OneDArray tot1(std::begin(t1), std::end(t1));

    xs->add(0, XS_t::TOTAL, tot1);

    double s1[][ng] = {{0.5}};

    XS_t::TwoDArray sct1(ng, ng, 0.0);

    for (int g = 0; g < ng; ++g)
    {
        for (int gp = 0; gp < ng; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
        }
    }

    xs->add(0, 0, sct1);

    double c[] = {1.0};
    double f[] = {0.2};
    double n[] = {0.4};

    XS_t::OneDArray chi(std::begin(c),std::end(c));
    XS_t::OneDArray fis(std::begin(f),std::end(f));
    XS_t::OneDArray nuf(std::begin(n),std::end(n));

    xs->add(0, XS_t::CHI,       chi);
    xs->add(0, XS_t::SIG_F,     fis);
    xs->add(0, XS_t::NU_SIG_F,  nuf);

    xs->complete();
    return xs;
}

Teuchos::RCP<XS_t> build_2grp_xs()
{
    constexpr int ng = 2;
    Teuchos::RCP<XS_t> xs = Teuchos::rcp(new XS_t());
    xs->set(0, ng);

    // make group boundaries
    XS_t::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 1.0;
    nbnd[2] = 0.00001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.0, 0.1};

    XS_t::OneDArray tot1(std::begin(t1), std::end(t1));

    xs->add(0, XS_t::TOTAL, tot1);

    double s1[][ng] = {{0.5, 0.0},{0.1, 0.05}};

    XS_t::TwoDArray sct1(ng, ng, 0.0);

    for (int g = 0; g < ng; ++g)
    {
        for (int gp = 0; gp < ng; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
        }
    }

    xs->add(0, 0, sct1);

    double c[] = {0.9, 0.10};
    double f[] = {0.2, 0.02};
    double n[] = {0.4, 0.05};

    XS_t::OneDArray chi(std::begin(c),std::end(c));
    XS_t::OneDArray fis(std::begin(f),std::end(f));
    XS_t::OneDArray nuf(std::begin(n),std::end(n));

    xs->add(0, XS_t::CHI,       chi);
    xs->add(0, XS_t::SIG_F,     fis);
    xs->add(0, XS_t::NU_SIG_F,  nuf);

    xs->complete();

    return xs;
}

}

void Keff_Tally_Tester::test_tally(int num_groups)
{
    int num_particles = 64;

    REQUIRE( num_groups==1 || num_groups==2 );

    Teuchos::RCP<XS_t> xs;
    if (num_groups == 1)
        xs = build_1grp_xs();
    else
        xs = build_2grp_xs();

    // Build geometry
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};
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
    sp_tally->begin_cycle(num_particles);
    tally.update_device();

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

    tally.update_host();
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
