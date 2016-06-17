//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/KCode_Solver_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  KCode_Solver_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KCode_Solver_Tester.hh"
#include "Test_XS.hh"
#include "../KCode_Solver.cuh"
#include "../Source_Transporter.hh"
#include "../Fission_Source.cuh"
#include "../Physics.cuh"
#include "../Cell_Tally.cuh"
#include "../Tallier.cuh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry       Geom;
typedef cuda_mc::Fission_Source<Geom>      Fission_Src;
typedef cuda_mc::Source_Transporter<Geom>  Transporter;
typedef cuda_mc::KCode_Solver<Geom>        Keff_Solver;

void KCode_Solver_Tester::test_transport(int num_groups)
{
    REQUIRE(num_groups==1 || num_groups==3 || num_groups==5);

    auto xs = Test_XS::build_xs(num_groups);

    std::vector<double> edges;
    std::vector<int> matids;

    if (num_groups==1)
    {
        edges = {0.0, 5.0, 10.0};
        matids = {0, 0, 0, 0, 0, 0, 0, 0};
    }
    else
    {
        edges = {0.0, 0.5, 1.0};
        matids = {0, 1, 1, 0, 0, 1, 1, 0};
    }

    def::size_type num_particles = 10000;

    // Build geometry
    auto geom = std::make_shared<Geom>(edges,edges,edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    pl->set("num_groups",xs->num_groups());
    pl->set("Np",num_particles);
    pl->set("implicit_capture",true);
    pl->set("variance reduction",std::string("roulette"));
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);
    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    cuda::Shared_Device_Ptr<Physics<Geom> > sdp_phys(phys);

    // Build cell tally
    std::cout << "Building Cell_Tally" << std::endl;
    auto sp_cell_tally = std::make_shared<Cell_Tally<Geom>>(
        sdp_geom,sdp_phys);
    std::vector<int> cells = {0, 1, 2, 3, 4, 5, 6, 7};
    sp_cell_tally->set_cells(cells);
    cuda::Shared_Device_Ptr<Cell_Tally<Geom> > cell_tally(sp_cell_tally);

    std::cout << "Building Tallier" << std::endl;
    auto tallier = std::make_shared<Tallier<Geom> >();
    tallier->add_cell_tally(cell_tally);

    // Build source
    auto source = std::make_shared<Fission_Src>(pl,sdp_geom,sdp_phys);

    // Build source transporter
    auto trans = std::make_shared<Transporter>(pl,sdp_geom,sdp_phys);

    Keff_Solver solver(pl);
    solver.set(trans,source,tallier);
    solver.solve();
    double keff = solver.keff();

    auto tally = sp_cell_tally->results();
    std::cout << "Tally result: ";
    for (auto x : tally)
        std::cout << x << " ";
    std::cout << std::endl;

    std::cout << "Computed keff: " << keff << std::endl;

    double exp_keff;
    if (num_groups==1)
        exp_keff = 0.697712;
    else if (num_groups==3)
        exp_keff = 1.014262;
    else if (num_groups==5)
        exp_keff = 0.392492;

    double tol = 10.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_SOFTEQ(exp_keff, keff, tol);

    // Statistical check on symmetry
    EXPECT_EQ( tally.size(), 8 );
    if (num_groups==1)
    {
        double mean = 0.0;
        for( auto x : tally )
            mean += x;
        mean /= static_cast<double>(tally.size());

        std::vector<double> exp(8,mean);
        EXPECT_VEC_SOFTEQ( exp, tally, tol );
    }
    else
    {
        double mean0 = 0.0;
        double mean1 = 0.0;
        int count0 = 0;
        int count1 = 0;
        for( int cell = 0; cell < tally.size(); ++cell )
        {
            if( matids[cell] == 0 )
            {
                mean0 += tally[cell];
                count0++;
            }
            else if( matids[cell] == 1 )
            {
                mean1 += tally[cell];
                count1++;
            }
        }
        mean0 /= static_cast<double>(count0);
        mean1 /= static_cast<double>(count1);

        std::vector<double> exp(8);
        for( int cell = 0; cell < matids.size(); ++cell )
        {
            if( matids[cell] == 0 )
                exp[cell] = mean0;
            else if( matids[cell] == 1 )
                exp[cell] = mean1;
        }
        EXPECT_VEC_SOFTEQ( exp, tally, tol );
    }
}

//---------------------------------------------------------------------------//
//                 end of KCode_Solver_Tester.cc

//---------------------------------------------------------------------------//
