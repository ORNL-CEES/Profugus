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
#include "cuda_rtk/RTK_Geometry.cuh"
#include "geometry/RTK_Geometry.hh"
#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace cuda_mc;


namespace
{

template <class Geom_DMM>
void run_transport(std::shared_ptr<Geom_DMM> geom_dmm,
                   int num_groups,
                   def::size_type num_particles,
                   def::size_type batch_size,
                   double &keff,
                   std::vector<double> &tally_mean,
                   std::vector<double> &tally_std_dev)
{
    typedef typename Geom_DMM::Geometry_t Geom;
    typedef cuda_mc::Fission_Source_DMM<Geom>  Fission_Src;
    typedef cuda_mc::Source_Transporter<Geom>  Transporter;
    typedef cuda_mc::KCode_Solver<Geom>        Keff_Solver;

    REQUIRE(num_groups==1 || num_groups==3 || num_groups==5);
    auto xs = Test_XS::build_xs(num_groups);

    // Build geometry
    auto sdp_geom = cuda::shared_device_ptr<Geom>(geom_dmm->device_instance());

    // Build physics
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    pl->set("num_groups",xs->num_groups());
    pl->set("Np",num_particles);
    pl->set("batch_size",batch_size);
    pl->set("implicit_capture",true);
    pl->set("variance reduction",std::string("roulette"));
    auto sdp_mat = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);
    auto phys = std::make_shared<Physics<Geom> >(pl,xs,sdp_mat);
    phys->set_geometry(sdp_geom);
    cuda::Shared_Device_Ptr<Physics<Geom> > sdp_phys(phys);

    pl->set("verbosity",std::string("medium"));

    // Build cell tally
    std::cout << "Building Cell_Tally" << std::endl;
    auto cell_tally = std::make_shared<Cell_Tally_DMM<Geom>>(
        sdp_geom,sdp_phys);
    int num_cells = geom_dmm->volumes().size();
    std::vector<int> cells;
    for (int cell = 0; cell < num_cells; ++cell)
        cells.push_back(cell);
    cell_tally->set_cells(cells,geom_dmm->volumes());

    std::cout << "Building Tallier" << std::endl;
    auto tallier = std::make_shared<Tallier_DMM<Geom>>();
    tallier->add_cell_tally(cell_tally);

    // Build source
    auto source = std::make_shared<Fission_Src>(pl,sdp_geom,sdp_phys,
            geom_dmm->lower(),geom_dmm->upper());

    // Build source transporter
    auto trans = std::make_shared<Transporter>(pl,sdp_geom,sdp_phys);

    Keff_Solver solver(pl);
    solver.set(trans,source,tallier);
    solver.solve();
    keff = solver.keff();

    tally_mean = cell_tally->results();
    std::cout << "Tally result: ";
    for (auto x : tally_mean)
        std::cout << x << " ";
    std::cout << std::endl;

    tally_std_dev = cell_tally->std_dev();
    std::cout << "Tally std dev: ";
    for (auto x : tally_std_dev)
        std::cout << x << " ";
    std::cout << std::endl;

    std::cout << "Computed keff: " << keff << std::endl;
}
}

void KCode_Solver_Tester::test_mesh(int num_groups)
{
    typedef cuda_profugus::Mesh_Geometry_DMM Geom_DMM;

    // Build different mesh depending on group structure
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

    // Build geometry
    auto geom_dmm = std::make_shared<Geom_DMM>(edges,edges,edges);
    geom_dmm->set_matids(matids);

    def::size_type num_particles = 10000;
    def::size_type batch_size    = num_particles / 5;

    // Run transport
    double keff;
    std::vector<double> tally_mean, tally_std_dev;
    run_transport(geom_dmm,num_groups,num_particles,batch_size,
                  keff,tally_mean,tally_std_dev);

    double exp_keff;
    if (num_groups==1)
        exp_keff = 0.697712;
    else if (num_groups==3)
        exp_keff = 1.014262;
    else if (num_groups==5)
        exp_keff = 0.392492;

    double tol = 20.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_SOFTEQ(exp_keff, keff, tol);

    // Statistical check on symmetry
    EXPECT_EQ( tally_mean.size(), 8 );
    if (num_groups==1)
    {
        double mean = 0.0;
        for( auto x : tally_mean )
            mean += x;
        mean /= static_cast<double>(tally_mean.size());

        std::vector<double> exp(8,mean);
        EXPECT_VEC_SOFTEQ( exp, tally_mean, tol );
    }
    else
    {
        double mean0 = 0.0;
        double mean1 = 0.0;
        int count0 = 0;
        int count1 = 0;
        for( int cell = 0; cell < tally_mean.size(); ++cell )
        {
            if( matids[cell] == 0 )
            {
                mean0 += tally_mean[cell];
                count0++;
            }
            else if( matids[cell] == 1 )
            {
                mean1 += tally_mean[cell];
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
        EXPECT_VEC_SOFTEQ( exp, tally_mean, tol );
    }
}

void KCode_Solver_Tester::test_rtk()
{
    typedef profugus::Core       Geometry_t;
    typedef Geometry_t::SP_Array SP_Core;
    typedef Geometry_t::Array_t  Core_t;
    typedef Core_t::SP_Object    SP_Lattice;
    typedef Core_t::Object_t     Lattice_t;
    typedef Lattice_t::SP_Object SP_Pin_Cell;
    typedef Lattice_t::Object_t  Pin_Cell_t;
    typedef Core_t::Vec_Int      Vec_Int;

    // make pin cells
    SP_Pin_Cell p1(std::make_shared<Pin_Cell_t>(1, 0.54, 0, 1.26, 14.28));
    SP_Pin_Cell p2(std::make_shared<Pin_Cell_t>(0, 1.26, 14.28));

    // make lattice
    SP_Lattice lat(std::make_shared<Lattice_t>(3, 3, 1, 2));

    // assign pins
    lat->assign_object(p1, 0); // fuel pins
    lat->assign_object(p2, 1); // guide tube

    // arrange pin-cells in lattice
    lat->id(0, 0, 0) = 0; // fuel pin
    lat->id(1, 0, 0) = 0; // fuel pin
    lat->id(2, 0, 0) = 0; // fuel pin
    lat->id(0, 1, 0) = 0; // fuel pin
    lat->id(1, 1, 0) = 1; // guide tube
    lat->id(2, 1, 0) = 0; // fuel pin
    lat->id(0, 2, 0) = 0; // fuel pin
    lat->id(1, 2, 0) = 0; // fuel pin
    lat->id(2, 2, 0) = 0; // fuel pin

    // complete lattice
    lat->complete(0.0, 0.0, 0.0);

    // make core
    SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
    core->assign_object(lat, 0);
    core->set_reflecting(Vec_Int(6, 1));
    core->complete(0.0, 0.0, 0.0);

    auto host_geometry = std::make_shared<Geometry_t>(core);

    // Now make Cuda geometry
    typedef cuda_profugus::Core_DMM Geom_DMM;
    auto geom_dmm = std::make_shared<Geom_DMM>(*host_geometry);

    def::size_type num_particles = 10000;
    def::size_type batch_size    = num_particles / 5;
    int num_groups = 3;

    // Run transport
    double keff;
    std::vector<double> tally_mean, tally_std_dev;
    run_transport(geom_dmm,num_groups,num_particles,batch_size,
                  keff,tally_mean,tally_std_dev);

    double tol = 20.0 / std::sqrt( static_cast<double>(num_particles) );
    double exp_keff = 1.184986; 
    EXPECT_SOFTEQ(exp_keff, keff, tol);
}

//---------------------------------------------------------------------------//
//                 end of KCode_Solver_Tester.cc

//---------------------------------------------------------------------------//
