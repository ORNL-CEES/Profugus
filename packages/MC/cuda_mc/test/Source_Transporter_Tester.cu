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
using profugus::XS;

namespace
{

Teuchos::RCP<XS> build_5grp_xs()
{
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, 5);

    std::vector<double> bnd(6, 0.0);
    bnd[0] = 100.0;
    bnd[1] = 10.0;
    bnd[2] = 1.0;
    bnd[3] = 0.1;
    bnd[4] = 0.01;
    bnd[5] = 0.001;
    xs->set_bounds(bnd);

    typename XS::OneDArray total(5);
    typename XS::TwoDArray scat(5, 5);

    double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
    double c[5] = {0.3770, 0.4421, 0.1809, 0.0, 0.0};
    double n[5] = {2.4*f[0], 2.4*f[1], 2.4*f[2], 2.4*f[3], 2.4*f[4]};
    typename XS::OneDArray fission(std::begin(f), std::end(f));
    typename XS::OneDArray chi(std::begin(c), std::end(c));
    typename XS::OneDArray nus(std::begin(n), std::end(n));
    xs->add(1, XS::SIG_F, fission);
    xs->add(1, XS::NU_SIG_F, nus);
    xs->add(1, XS::CHI, chi);

    // mat 0
    total[0] = 5.2 ;
    total[1] = 11.4;
    total[2] = 18.2;
    total[3] = 29.9;
    total[4] = 27.3;
    xs->add(0, XS::TOTAL, total);

    // mat 1
    total[0] = 5.2  + f[0];
    total[1] = 11.4 + f[1];
    total[2] = 18.2 + f[2];
    total[3] = 29.9 + f[3];
    total[4] = 27.3 + f[4];
    xs->add(1, XS::TOTAL, total);

    scat(0, 0) = 1.2;
    scat(1, 0) = 0.9;
    scat(1, 1) = 3.2;
    scat(2, 0) = 0.4;
    scat(2, 1) = 2.8;
    scat(2, 2) = 6.9;
    scat(2, 3) = 1.5;
    scat(3, 0) = 0.1;
    scat(3, 1) = 2.1;
    scat(3, 2) = 5.5;
    scat(3, 3) = 9.7;
    scat(3, 4) = 2.1;
    scat(4, 1) = 0.2;
    scat(4, 2) = 1.3;
    scat(4, 3) = 6.6;
    scat(4, 4) = 9.9;
    xs->add(0, 0, scat);
    xs->add(1, 0, scat);

    xs->complete();
    return xs;
}

Teuchos::RCP<XS> build_3grp_xs()
{
    const int ng = 3;
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, ng);

    // make group boundaries
    XS::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 1.0;
    nbnd[2] = 0.01;
    nbnd[3] = 0.0001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.1, 1.6, 2.9};
    double t2[ng] = {10.0, 11.3, 16.2};

    XS::OneDArray tot1(std::begin(t1), std::end(t1));
    XS::OneDArray tot2(std::begin(t2), std::end(t2));

    xs->add(0, XS::TOTAL, tot1);
    xs->add(1, XS::TOTAL, tot2);

    double s1[][3] = {{0.7, 0.0, 0.0},
                      {0.2, 0.3, 0.0},
                      {0.1, 0.7, 1.9}};

    double s2[][3] = {{2.7, 0.0, 0.0},
                      {2.2, 2.3, 0.0},
                      {2.1, 2.7, 3.9}};

    XS::TwoDArray sct1(ng, ng, 0.0);
    XS::TwoDArray sct2(ng, ng, 0.0);

    for (int g = 0; g < 3; ++g)
    {
        for (int gp = 0; gp < 3; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
            sct2(g, gp) = s2[g][gp];
        }
    }

    xs->add(0, 0, sct1);
    xs->add(1, 0, sct2);

    double c2[] = {0.4, 0.6, 0.0};
    double f2[] = {3.2, 4.2, 0.0};
    double n2[] = {2.4*3.2, 2.4*4.2, 0.0};

    XS::OneDArray chi2(std::begin(c2), std::end(c2));
    XS::OneDArray fis2(std::begin(f2), std::end(f2));
    XS::OneDArray nuf2(std::begin(n2), std::end(n2));

    xs->add(1, XS::CHI, chi2);
    xs->add(1, XS::SIG_F, fis2);
    xs->add(1, XS::NU_SIG_F, nuf2);

    xs->complete();
    return xs;
}

Teuchos::RCP<XS> build_1grp1mat_xs()
{
    const int ng = 1;
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, ng);

    // make group boundaries
    XS::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 0.00001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.0};

    XS::OneDArray tot1(std::begin(t1), std::end(t1));

    xs->add(0, XS::TOTAL, tot1);

    double s1[][3] = {{0.9}};

    XS::TwoDArray sct1(ng, ng, 0.0);

    for (int g = 0; g < ng; ++g)
    {
        for (int gp = 0; gp < ng; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
        }
    }

    xs->add(0, 0, sct1);

    xs->complete();
    return xs;
}

}

void Source_Transporter_Tester::test_transport(int num_groups)
{
    REQUIRE(num_groups==1 || num_groups==3 || num_groups==5);

    Teuchos::RCP<XS> xs;
    if (num_groups==1)
        xs = build_1grp1mat_xs();
    else if (num_groups==3)
        xs = build_3grp_xs();
    else
        xs = build_5grp_xs();

    int Np = 10000;

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
