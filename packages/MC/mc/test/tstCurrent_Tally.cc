//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstCurrent_Tally.cc
 * \author Steven Hamilton
 * \date   Fri Aug 21 15:46:50 2015
 * \brief  Current_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <utility>
#include <algorithm>
#include <memory>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "geometry/Mesh_Geometry.hh"
#include "../Current_Tally.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class CurrentTallyTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Mesh_Geometry             Geometry_t;
    typedef profugus::Current_Tally<Geometry_t> Current_Tally;
    typedef std::shared_ptr<Geometry_t>         SP_Geometry;
    typedef profugus::Physics<Geometry_t>       Physics_t;
    typedef std::shared_ptr<Physics_t>          SP_Physics;
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

        d_db = Teuchos::rcp(new ParameterList_t("test"));
        d_db->set("problem_name",std::string("current_tally_test"));

        d_x_edges = {0.0, 0.5, 1.0};
        d_y_edges = {2.0, 3.0, 4.0, 5.0};
        d_z_edges = {-0.5, 0.0, 0.5};

        build_geometry();
        build_physics();

        d_tally = std::make_shared<Current_Tally>(
            d_db,d_physics,d_x_edges,d_y_edges,d_z_edges);
    }

    void build_geometry()
    {
        d_geometry = std::make_shared<Geometry_t>(
            d_x_edges,d_y_edges,d_z_edges);
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

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));

        xs->add(0, XS_t::TOTAL, tot1);

        double s1[][3] = {{0.7, 0.0, 0.0},
                          {0.2, 0.3, 0.0},
                          {0.1, 0.7, 1.9}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);

        for (int g = 0; g < 3; ++g)
        {
            for (int gp = 0; gp < 3; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
            }
        }

        xs->add(0, 0, sct1);

        xs->complete();

        d_physics = std::make_shared<Physics_t>(d_db, xs);
        d_physics->set_geometry(d_geometry);
    }

  protected:
    // >>> DATA

    RCP_Std_DB  d_db;
    SP_Geometry d_geometry;
    SP_Physics  d_physics;

    std::vector<double> d_x_edges;
    std::vector<double> d_y_edges;
    std::vector<double> d_z_edges;

    std::shared_ptr<Current_Tally> d_tally;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CurrentTallyTest, multi_cell)
{
    // d_x_edges = {0.0, 0.5, 1.0};
    // d_y_edges = {2.0, 3.0, 4.0, 5.0};
    // d_z_edges = {-0.5, 0.0, 0.5};

    Particle_t p;
    p.set_wt(1.0);

    //
    // History 1
    //

    // Tally on x surface 1
    d_geometry->initialize({0.5, 2.5, -0.2}, {1.0, 0.0, 0.0}, p.geo_state());
    d_tally->tally_surface(p);

    // Again on same surface
    d_geometry->initialize({0.5, 2.2, -0.4}, {-0.5, 0.0, 0.5}, p.geo_state());
    d_tally->tally_surface(p);

    // Not on surface, should be no contribution
    d_geometry->initialize({0.6, 3.4, 0.4}, {0.2, 0.5, -0.1}, p.geo_state());
    d_tally->tally_surface(p);

    d_tally->end_history();

    //
    // History 2
    //

    p.set_wt(0.5);

    // Tally on z surface 5
    d_geometry->initialize({0.7, 4.4, -0.5}, {0.0, 1.0, -1.0}, p.geo_state());
    d_tally->tally_surface(p);

    // Not on surface, should be no contribution
    d_geometry->initialize({0.1, 4.1, -0.2}, {0.1, -0.2, 0.8}, p.geo_state());
    d_tally->tally_surface(p);

    d_tally->end_history();

    //
    // History 3
    //

    p.set_wt(2.0);

    // Tally on same x surface as particle 1
    d_geometry->initialize({0.5, 2.1, -0.5}, {1.0, 0.0, 0.0}, p.geo_state());
    d_tally->tally_surface(p);

    // Tally on y surface 2
    d_geometry->initialize({0.4, 3.0, -0.1}, {1.0, 1.0, 0.0}, p.geo_state());
    d_tally->tally_surface(p);

    d_tally->end_history();

    d_tally->finalize(3.0 * static_cast<double>(nodes));

    // Get the current results
    auto x_current = d_tally->x_current();
    auto y_current = d_tally->y_current();
    auto z_current = d_tally->z_current();
    auto x_current_std_dev = d_tally->x_current_std_dev();
    auto y_current_std_dev = d_tally->y_current_std_dev();
    auto z_current_std_dev = d_tally->z_current_std_dev();

    EXPECT_EQ( 18, x_current.size() );
    EXPECT_EQ( 16, y_current.size() );
    EXPECT_EQ( 18, z_current.size() );
    EXPECT_EQ( 18, x_current_std_dev.size() );
    EXPECT_EQ( 16, y_current_std_dev.size() );
    EXPECT_EQ( 18, z_current_std_dev.size() );

    double x_area = 0.5;
    double y_area = 0.5 * 0.5;
    double z_area = 0.5;

    // Check value on x face 1
    double expected =  1.0 - 1.0  + // History 1
                       2.0;         // History 3
    expected /= (3.0  * x_area);
    EXPECT_SOFT_EQ(expected, x_current[1]);

    // Check value on y face
    expected = 2.0; // History 3
    expected /= (3.0 * y_area);
    EXPECT_SOFT_EQ(expected, y_current[2]);

    // Check value on z face
    expected = -0.5; // History 2
    expected /= (3.0 * z_area);
    EXPECT_SOFT_EQ(expected, z_current[5]);

    if (nodes == 1)
    {
        expected = 2.309401076758503 / std::sqrt(3.0);
        EXPECT_SOFT_EQ(expected, x_current_std_dev[1]);
        expected = 4.618802153517007 / std::sqrt(3.0);
        EXPECT_SOFT_EQ(expected, y_current_std_dev[2]);
        expected = 0.577350269189626 / std::sqrt(3.0);
        EXPECT_SOFT_EQ(expected, z_current_std_dev[5]);
    }
    else if (nodes == 4)
    {
        expected = 1.969463855669324 / std::sqrt(12.0);
        EXPECT_SOFT_EQ(expected, x_current_std_dev[1]);
        expected = 3.938927711338648 / std::sqrt(12.0);
        EXPECT_SOFT_EQ(expected, y_current_std_dev[2]);
        expected = 0.492365963917331 / std::sqrt(12.0);
        EXPECT_SOFT_EQ(expected, z_current_std_dev[5]);
    }

    if (node == 0)
    {
        std::cout << "X currents: ";
        for( auto &x : x_current )
            std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Y currents: ";
        for( auto &x : y_current )
            std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Z currents: ";
        for( auto &x : z_current )
            std::cout << x << " ";
        std::cout << std::endl;
    }

    // Get the flux results
    auto x_flux = d_tally->x_flux();
    auto y_flux = d_tally->y_flux();
    auto z_flux = d_tally->z_flux();
    auto x_flux_std_dev = d_tally->x_flux_std_dev();
    auto y_flux_std_dev = d_tally->y_flux_std_dev();
    auto z_flux_std_dev = d_tally->z_flux_std_dev();

    EXPECT_EQ( 18, x_flux.size() );
    EXPECT_EQ( 16, y_flux.size() );
    EXPECT_EQ( 18, z_flux.size() );
    EXPECT_EQ( 18, x_flux_std_dev.size() );
    EXPECT_EQ( 16, y_flux_std_dev.size() );
    EXPECT_EQ( 18, z_flux_std_dev.size() );

    // Check value on x face 1
    expected =  1.0 + 1.0*std::sqrt(2.0)  + // History 1
                      2.0;   // History 3
    expected /= (3.0 * x_area);
    EXPECT_SOFT_EQ(expected, x_flux[1]);

    // Check value on y face
    expected = 2.0*std::sqrt(2.0); // History 3
    expected /= (3.0 * y_area);
    EXPECT_SOFT_EQ(expected, y_flux[2]);

    // Check value on z face
    expected = 0.5*std::sqrt(2.0); // History 2
    expected /= (3.0 * z_area);
    EXPECT_SOFT_EQ(expected, z_flux[5]);

    if (nodes == 1)
    {
        expected = 2.581988897471611 / std::sqrt(3.0);
        EXPECT_SOFT_EQ(expected, x_flux_std_dev[1]);
        expected = 6.531972647421808 / std::sqrt(3.0);
        EXPECT_SOFT_EQ(expected, y_flux_std_dev[2]);
        expected = 0.816496580927726 / std::sqrt(3.0);
        EXPECT_SOFT_EQ(expected, z_flux_std_dev[5]);
    }
    else if (nodes == 4)
    {
        expected = 2.201927530252721 / std::sqrt(12.0);
        EXPECT_SOFT_EQ(expected, x_flux_std_dev[1]);
        expected = 5.570484990582331 / std::sqrt(12.0);
        EXPECT_SOFT_EQ(expected, y_flux_std_dev[2]);
        expected = 0.696310623822791 / std::sqrt(12.0);
        EXPECT_SOFT_EQ(expected, z_flux_std_dev[5]);
    }

    if (node == 0)
    {
        std::cout << "X fluxes: ";
        for( auto &x : x_flux )
            std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Y fluxes: ";
        for( auto &x : y_flux )
            std::cout << x << " ";
        std::cout << std::endl;
        std::cout << "Z fluxes: ";
        for( auto &x : z_flux )
            std::cout << x << " ";
        std::cout << std::endl;
    }
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstCurrent_Tally.cc
//
//---------------------------------------------------------------------------//
