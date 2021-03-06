//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstFission_Tally.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:46:50 2015
 * \brief  Fission_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <utility>
#include <algorithm>
#include <memory>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "geometry/RTK_Geometry.hh"
#include "../Fission_Tally.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class FissionTallyTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Core                      Geometry_t;
    typedef profugus::Fission_Tally<Geometry_t> Fission_Tally;
    typedef Fission_Tally::SP_Geometry          SP_Geometry;
    typedef Fission_Tally::Physics_t            Physics_t;
    typedef Fission_Tally::SP_Physics           SP_Physics;
    typedef Physics_t::Particle_t               Particle_t;
    typedef Physics_t::XS_t                     XS_t;
    typedef Physics_t::RCP_XS                   RCP_XS;

    typedef Teuchos::ParameterList        ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t> RCP_Std_DB;

  protected:
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();

        db = Teuchos::rcp(new ParameterList_t("test"));
        db->set("problem_name",std::string("fission_tally_test"));

        build_geometry();
        build_physics();

        tally = std::make_shared<Fission_Tally>(physics);

        // Make Cartesian mesh
        std::vector<double> xy_edges = {0.0, 10.0, 20.0};
        std::vector<double> z_edges  = {0.0, 20.0};
        auto mesh = std::make_shared<profugus::Cartesian_Mesh>(
            xy_edges,xy_edges,z_edges);
        tally->set_mesh(mesh);
    }

    void build_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // make boxes (10x10x20)
        SP_Pin_Cell p0(std::make_shared<Pin_Cell_t>(0, 10.0, 20.0));
        SP_Pin_Cell p1(std::make_shared<Pin_Cell_t>(1, 10.0, 20.0));

        // make lattice
        /*
          1 0
          0 1
         */
        SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 2));

        // assign pins
        lat->assign_object(p0, 0); // mat 0 box
        lat->assign_object(p1, 1); // mat 1 box

        // arrange pin-cells in lattice
        lat->id(0, 0, 0) = 0; // mat 0 box
        lat->id(1, 0, 0) = 1; // mat 1 box
        lat->id(0, 1, 0) = 1; // mat 1 box
        lat->id(1, 1, 0) = 0; // mat 0 box

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->complete(0.0, 0.0, 0.0);

        geometry = std::make_shared<Geometry_t>(core);
        EXPECT_EQ(4, core->num_cells());
    }

    void build_physics()
    {
        const int ng = 3;

        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(4, 0.0);
        nbnd[0] = 100.0; nbnd[1] = 1.0; nbnd[2] = 0.01; nbnd[3] = 0.0001;
        xs->set_bounds(nbnd);

        double t1[] = {1.1, 1.6, 2.9};
        double t2[] = {10.0, 11.3, 16.2};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));
        XS_t::OneDArray tot2(std::begin(t2), std::end(t2));

        xs->add(0, XS_t::TOTAL, tot1);
        xs->add(1, XS_t::TOTAL, tot2);

        double s1[][3] = {{0.7, 0.0, 0.0},
                          {0.2, 0.3, 0.0},
                          {0.1, 0.7, 1.9}};

        double s2[][3] = {{2.7, 0.0, 0.0},
                          {2.2, 2.3, 0.0},
                          {2.1, 2.7, 3.9}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);
        XS_t::TwoDArray sct2(ng, ng, 0.0);

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

        XS_t::OneDArray chi2(std::begin(c2), std::end(c2));
        XS_t::OneDArray fis2(std::begin(f2), std::end(f2));
        XS_t::OneDArray nuf2(std::begin(n2), std::end(n2));

        xs->add(1, XS_t::CHI, chi2);
        xs->add(1, XS_t::SIG_F, fis2);
        xs->add(1, XS_t::NU_SIG_F, nuf2);

        xs->complete();

        physics = std::make_shared<Physics_t>(db, xs);
        physics->set_geometry(geometry);
    }

  protected:
    // >>> DATA

    RCP_Std_DB  db;
    SP_Geometry geometry;
    SP_Physics  physics;

    std::shared_ptr<Fission_Tally> tally;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(FissionTallyTest, multi_cell)
{
    Particle_t p;
    p.set_wt(1.0);
    p.set_group(0);

    // History 1

    geometry->initialize({1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(0, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // tally here -- non-fissionable cell, no contribution
    tally->accumulate(3.0, p);

    geometry->initialize({11.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(1, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // tally here
    tally->accumulate(3.0, p);

    geometry->initialize({5.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(2, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // no tally here
    tally->accumulate(3.0, p);

    geometry->initialize({15.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(3, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // tally here -- non-fissionable cell, no contribution
    tally->accumulate(4.0, p);
    tally->accumulate(2.0, p);

    geometry->initialize({9.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(2, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // no tally here
    tally->accumulate(3.0, p);

    tally->end_history();

    // History 2

    p.set_group(1);

    geometry->initialize({5.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(2, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // no tally here
    tally->accumulate(3.0, p);
    tally->accumulate(6.0, p);
    tally->accumulate(9.0, p);

    geometry->initialize({15.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(3, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // tally here -- non-fissionable cell, no contribution
    tally->accumulate(2.0, p);
    tally->accumulate(1.0, p);

    tally->end_history();

    // History 3

    p.set_group(0);

    geometry->initialize({5.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(2, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // no tally here
    tally->accumulate(7.0, p);
    tally->accumulate(6.0, p);

    geometry->initialize({15.0, 15.0, 1.0}, {1.0, 1.0, 1.0}, p.geo_state());
    EXPECT_EQ(3, geometry->cell(p.geo_state()));
    p.set_matid(geometry->matid(p.geo_state()));

    // tally here
    tally->accumulate(8.0, p);

    tally->end_history();

    // Finalize
    tally->finalize(3 * nodes);

    // Get the results
    auto results = tally->results();
    double tol = 1.0e-8;

    EXPECT_EQ(4, results.size());

    const auto &r0 = results[0];
    const auto &r1 = results[1];
    const auto &r2 = results[2];
    const auto &r3 = results[3];

    // Cells 0 and 3 have no fissionable material
    EXPECT_SOFTEQ( 0.0, r0.first,  tol );
    EXPECT_SOFTEQ( 0.0, r0.second, tol );
    EXPECT_SOFTEQ( 0.0, r3.first,  tol );
    EXPECT_SOFTEQ( 0.0, r3.second, tol );

    EXPECT_SOFTEQ(3.840e-3,  r1.first, tol);
    EXPECT_SOFTEQ(5.456e-2,  r2.first, tol);

    if( nodes == 1 )
    {
        EXPECT_SOFTEQ(3.840e-3,          r1.second, tol);
        EXPECT_SOFTEQ(0.019674796059934, r2.second, tol);
    }
    if( nodes == 4 )
    {
        EXPECT_SOFTEQ(0.001637381501611, r1.second, tol);
        EXPECT_SOFTEQ(0.008389361228257, r2.second, tol);
    }
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstFission_Tally.cc
//---------------------------------------------------------------------------//
