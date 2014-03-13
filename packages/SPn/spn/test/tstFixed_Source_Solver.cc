//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstFixed_Source_Solver.cc
 * \author Thomas M. Evans
 * \date   Sun Nov 18 19:05:06 2012
 * \brief  Test of Fixed_Source_Solver class.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "gtest/utils_gtest.hh"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "utils/Definitions.hh"
#include "mesh/Partitioner.hh"
#include "../Dimensions.hh"
#include "../Fixed_Source_Solver.hh"
#include "Test_XS.hh"

using namespace std;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Inf_Med_Solver_FVTest : public testing::Test
{
  protected:

    typedef profugus::Fixed_Source_Solver          Fixed_Source_Solver;
    typedef Teuchos::RCP<Fixed_Source_Solver>      RCP_Linear_Solver;
    typedef Fixed_Source_Solver::RCP_ParameterList RCP_ParameterList;
    typedef Fixed_Source_Solver::RCP_Mat_DB        RCP_Mat_DB;
    typedef Fixed_Source_Solver::RCP_Dimensions    RCP_Dimensions;
    typedef Fixed_Source_Solver::RCP_Mesh          RCP_Mesh;
    typedef Fixed_Source_Solver::RCP_Indexer       RCP_Indexer;
    typedef Fixed_Source_Solver::RCP_Global_Data   RCP_Global_Data;
    typedef Fixed_Source_Solver::Matrix_t          Matrix_t;
    typedef Fixed_Source_Solver::Vector_t          Vector_t;
    typedef Fixed_Source_Solver::Linear_System_t   Linear_System_t;
    typedef profugus::Partitioner                  Partitioner;
    typedef Fixed_Source_Solver::State_t           State;
    typedef Teuchos::RCP<State>                    RCP_State;
    typedef Linear_System_t::Array_Dbl             Array_Dbl;
    typedef Fixed_Source_Solver::External_Source   External_Source;
    typedef External_Source::Source_Shapes         Source_Shapes;
    typedef External_Source::Shape                 Shape;
    typedef External_Source::Source_Field          Source_Field;
    typedef External_Source::ID_Field              ID_Field;

  protected:

    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

    void build(int order,
               int Ng)
    {
        num_groups = Ng;
        eqn_order  = order;

        // build 4x4x4 mesh
        RCP_ParameterList db = Teuchos::rcp(new Teuchos::ParameterList("test"));

        if (cx.size() == 0)
        {
            db->set("delta_x", 1.0);
            db->set("delta_y", 1.0);
            db->set("delta_z", 1.0);

            db->set("num_cells_i", 3);
            db->set("num_cells_j", 3);
            db->set("num_cells_k", 3);
        }
        else
        {
            db->set("x_edges", cx);
            db->set("y_edges", cy);
            db->set("z_edges", cz);
        }

        if (nodes == 2)
        {
            db->set("num_blocks_i", 2);
        }
        if (nodes == 4)
        {
            db->set("num_blocks_i", 2);
            db->set("num_blocks_j", 2);
        }

        db->set("tolerance",1.0e-8);

        Partitioner p(db);
        p.build();

        mesh    = p.get_mesh();
        indexer = p.get_indexer();
        data    = p.get_global_data();

        solver  = Teuchos::rcp(new Fixed_Source_Solver(db));

        if (num_groups == 1)
            mat = one_grp::make_mat(3, mesh->num_cells());
        else if (num_groups == 2)
            mat = two_grp::make_mat(3, mesh->num_cells());
        else
            mat = three_grp::make_mat(3, mesh->num_cells());

        EXPECT_FALSE(mat.is_null());
        EXPECT_FALSE(mesh.is_null());
        EXPECT_FALSE(indexer.is_null());
        EXPECT_FALSE(data.is_null());

        bool success = true, verbose = true;
        try
        {
            dim = Teuchos::rcp(new profugus::Dimensions(eqn_order));

            solver->setup(dim, mat, mesh, indexer, data);
        }
        TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

        // make a state object
        state = Teuchos::rcp(new State(mesh, Ng));
    }

  protected:

    RCP_Mesh        mesh;
    RCP_Indexer     indexer;
    RCP_Global_Data data;

    RCP_Mat_DB     mat;
    RCP_Dimensions dim;

    RCP_Linear_Solver solver;

    RCP_State state;

    int num_groups, eqn_order;

    int node, nodes;

    Array_Dbl cx, cy, cz;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Inf_Med_Solver_FVTest, 1Grp_SP1)
{
    build(1, 1);
    EXPECT_EQ(1, dim->num_equations());

    // make the source
    External_Source q(mesh->num_cells());
    {
        Source_Shapes shapes(1, Shape(1, 1.2));
        ID_Field srcids(mesh->num_cells(), 0);
        Source_Field source(mesh->num_cells(), 1.0);

        q.set(srcids, shapes, source);
    }

    solver->solve(q);

    const Vector_t &x = solver->get_LHS();
    for (int i = 0; i < x.MyLength(); ++i)
    {
        EXPECT_SOFTEQ(2.0, x[i], 1.0e-6);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Inf_Med_Solver_FVTest, 1Grp_SP3)
{
    build(3, 1);
    EXPECT_EQ(2, dim->num_equations());

    // make the source
    External_Source q(mesh->num_cells());
    {
        Source_Shapes shapes(1, Shape(1, 1.2));
        ID_Field srcids(mesh->num_cells(), 0);
        Source_Field source(mesh->num_cells(), 1.0);

        q.set(srcids, shapes, source);
    }

    solver->solve(q);

    const Vector_t &x = solver->get_LHS();
    for (int cell = 0; cell < mesh->num_cells(); ++cell)
    {
        double phi_0 = x[cell*2] - 2.0/3.0 * x[1 + cell*2];
        EXPECT_SOFTEQ(2.0, phi_0, 1.0e-6);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Inf_Med_Solver_FVTest, 3Grp_SP1)
{
    build(1, 3);
    EXPECT_EQ(1, dim->num_equations());

    // make the source
    External_Source q(mesh->num_cells());
    {
        Source_Shapes shapes(1, Shape(3, 0.0));
        shapes[0][0] = 1.2;
        shapes[0][1] = 1.3;
        shapes[0][2] = 1.4;

        ID_Field srcids(mesh->num_cells(), 0);
        Source_Field source(mesh->num_cells(), 1.0);

        q.set(srcids, shapes, source);
    }

    solver->solve(q);

    const Linear_System_t &system = solver->get_linear_system();

    const Vector_t &x = solver->get_LHS();
    for (int cell = 0; cell < mesh->num_cells(); ++cell)
    {
        int g0       = 0 + 0 * 3 + cell * 3 * 1;
        int g1       = 1 + 0 * 3 + cell * 3 * 1;
        int g2       = 2 + 0 * 3 + cell * 3 * 1;
        EXPECT_EQ(g0, system.index(0, 0, cell));
        EXPECT_EQ(g1, system.index(1, 0, cell));
        EXPECT_EQ(g2, system.index(2, 0, cell));

        double phi_0 = x[g0];
        double phi_1 = x[g1];
        double phi_2 = x[g2];
        EXPECT_SOFTEQ(23.376775173864782, phi_0, 1.0e-6);
        EXPECT_SOFTEQ(26.285032257831212, phi_1, 1.0e-6);
        EXPECT_SOFTEQ(21.044148232485092, phi_2, 1.0e-6);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Inf_Med_Solver_FVTest, 3Grp_SP5)
{
    build(5, 3);
    EXPECT_EQ(3, dim->num_equations());

    // make the source
    External_Source q(mesh->num_cells());
    {
        Source_Shapes shapes(1, Shape(3, 0.0));
        shapes[0][0] = 1.2;
        shapes[0][1] = 1.3;
        shapes[0][2] = 1.4;

        ID_Field srcids(mesh->num_cells(), 0);
        Source_Field source(mesh->num_cells(), 1.0);

        q.set(srcids, shapes, source);
    }

    solver->solve(q);

    const Linear_System_t &system = solver->get_linear_system();

    const Vector_t &x = solver->get_LHS();
    double eps = 1.0e-4;
    for (int cell = 0; cell < mesh->num_cells(); ++cell)
    {
        int g00 = 0 + 0 * 3 + cell * 3 * 3;
        int g10 = 1 + 0 * 3 + cell * 3 * 3;
        int g20 = 2 + 0 * 3 + cell * 3 * 3;
        EXPECT_EQ(g00, system.index(0, 0, cell));
        EXPECT_EQ(g10, system.index(1, 0, cell));
        EXPECT_EQ(g20, system.index(2, 0, cell));

        int g01 = 0 + 1 * 3 + cell * 3 * 3;
        int g11 = 1 + 1 * 3 + cell * 3 * 3;
        int g21 = 2 + 1 * 3 + cell * 3 * 3;
        EXPECT_EQ(g01, system.index(0, 1, cell));
        EXPECT_EQ(g11, system.index(1, 1, cell));
        EXPECT_EQ(g21, system.index(2, 1, cell));

        int g02 = 0 + 2 * 3 + cell * 3 * 3;
        int g12 = 1 + 2 * 3 + cell * 3 * 3;
        int g22 = 2 + 2 * 3 + cell * 3 * 3;
        EXPECT_EQ(g02, system.index(0, 2, cell));
        EXPECT_EQ(g12, system.index(1, 2, cell));
        EXPECT_EQ(g22, system.index(2, 2, cell));

        double phi_0 = x[g00] - 2.0/3.0*x[g01] + 8.0/15.0*x[g02];
        double phi_1 = x[g10] - 2.0/3.0*x[g11] + 8.0/15.0*x[g12];
        double phi_2 = x[g20] - 2.0/3.0*x[g21] + 8.0/15.0*x[g22];
        EXPECT_SOFTEQ(23.376775173864782, phi_0, eps);
        EXPECT_SOFTEQ(26.285032257831212, phi_1, eps);
        EXPECT_SOFTEQ(21.044148232485092, phi_2, eps);
    }

    // fill the state
    solver->write_state(*state);
    State::View_Field phi = state->flux();
    EXPECT_EQ(mesh->num_cells() * 3, phi.size());

    for (int cell = 0; cell < mesh->num_cells(); ++cell)
    {
        int g0 = cell + 0 * mesh->num_cells();
        int g1 = cell + 1 * mesh->num_cells();
        int g2 = cell + 2 * mesh->num_cells();

        EXPECT_SOFTEQ(23.376775173864782, phi[g0], eps);
        EXPECT_SOFTEQ(26.285032257831212, phi[g1], eps);
        EXPECT_SOFTEQ(21.044148232485092, phi[g2], eps);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstFixed_Source_Solver.cc
//---------------------------------------------------------------------------//
