//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source_Transporter_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Source_Transporter_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source_Transporter_Tester.hh"
#include "Test_XS.hh"
#include "../Source_Transporter.hh"
#include "../Uniform_Source.cuh"
#include "../Physics.cuh"
#include "../Cell_Tally.cuh"
#include "../Tallier.cuh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;

typedef cuda_profugus::Mesh_Geometry      Geom;
typedef cuda_mc::Uniform_Source<Geom>     Uniform_Src;
typedef cuda_mc::Source_Transporter<Geom> Transporter;

void Source_Transporter_Tester::test_transport(int num_groups)
{
    auto xs = Test_XS::build_xs(num_groups);

    def::size_type Np = 10000;

    // Build geometry
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids;
    if (num_groups == 1)
        matids = {0, 0, 0, 0, 0, 0, 0, 0};
    else
        matids = {0, 1, 1, 0, 0, 1, 1, 0};
    auto geom = std::make_shared<Geom>(edges,edges,edges);
    geom->set_matids(matids);
    cuda::Shared_Device_Ptr<cuda_profugus::Mesh_Geometry> sdp_geom(geom);

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    pl->set("num_groups",xs->num_groups());
    pl->set("Np",Np);
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
    auto sp_tallier = std::make_shared<Tallier<Geom> >();
    sp_tallier->add_cell_tally(cell_tally);
    cuda::Shared_Device_Ptr<Tallier<Geom>> tallier(sp_tallier);

    // Build box shape for source
    std::vector<double> src_bounds = {edges.front(), edges.back(),
                                      edges.front(), edges.back(),
                                      edges.front(), edges.back()};
    REQUIRE( src_bounds.size() == 6 );
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            src_bounds[0], src_bounds[1],
            src_bounds[2], src_bounds[3],
            src_bounds[4], src_bounds[5]);

    // Build source
    auto source = std::make_shared<Uniform_Src>(pl,sdp_geom);
    source->build_source(src_shape);

    // Build source transporter
    Transporter trans(pl,sdp_geom,sdp_phys);
    trans.set(tallier);
    trans.solve(source);

    sp_tallier->finalize(Np);
    auto tally = sp_cell_tally->results();
    std::cout << "Tally result: ";
    for( auto x : tally )
        std::cout << x << " ";
    std::cout << std::endl;

    // Test statistics on output using symmetry
    EXPECT_EQ( tally.size(), 8 );

    if (num_groups==1)
    {
        double mean = 0.0;
        for( auto x : tally )
            mean += x;
        mean /= static_cast<double>(tally.size());

        double tol = 10.0 / std::sqrt( static_cast<double>(Np) );

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

        double tol = 10.0 / std::sqrt( static_cast<double>(Np) );

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
//                 end of Source_Transporter_Tester.cc

//---------------------------------------------------------------------------//
