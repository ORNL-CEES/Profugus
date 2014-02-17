//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstState.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 17 14:01:48 2014
 * \brief  State unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Teuchos_RCP.hpp"

#include "mesh/Partitioner.hh"
#include "../State.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class State_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::State                State;
    typedef Teuchos::RCP<State>            RCP_State;
    typedef profugus::Partitioner          Partitioner;
    typedef State::const_View_Field        const_View_Field;
    typedef State::View_Field              View_Field;
    typedef Partitioner::RCP_Mesh          RCP_Mesh;
    typedef Partitioner::RCP_Indexer       RCP_Indexer;
    typedef Partitioner::RCP_Global_Data   RCP_Global_Data;
    typedef Partitioner::RCP_ParameterList RCP_ParameterList;
    typedef Partitioner::ParameterList     ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        pl = Teuchos::rcp(new ParameterList("Part"));
    }

    void init(int Pi, int Pj, int Ng)
    {
        pl->set("num_blocks_i", Pi);
        pl->set("num_blocks_j", Pj);
        pl->set("num_cells_i", 5);
        pl->set("num_cells_j", 8);
        pl->set("num_cells_k", 3);
        pl->set("delta_x", 0.1);
        pl->set("delta_y", 0.1);
        pl->set("delta_z", 0.1);
        pl->set("num_z_blocks", 1);

        Partitioner p(pl);
        p.build();

        mesh    = p.get_mesh();
        indexer = p.get_indexer();
        gdata   = p.get_global_data();
        state   = Teuchos::rcp(new State(mesh, Ng));
    }

  protected:
    // >>> Data that get re-initialized between tests

    int node, nodes;

    RCP_State state;

    RCP_ParameterList pl;
    RCP_Mesh          mesh;
    RCP_Indexer       indexer;
    RCP_Global_Data   gdata;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(State_Test, 1PE_1Grp)
{
    if (nodes != 1)
        return;

    init(1, 1, 1);

    EXPECT_EQ(120, state->size());
    EXPECT_EQ(120, state->num_cells());
    EXPECT_EQ(1, state->num_groups());

    View_Field md = state->flux();
    for (int c = 0; c < 120; ++c)
    {
        EXPECT_EQ(0.0, md[c]);
    }

    const State &const_state = *state;
    const_View_Field cd      = const_state.flux();
    for (int c = 0; c < 120; ++c)
    {
        EXPECT_EQ(0.0, cd[c]);
    }

    md[1] = 1.0;
    md[5] = 2.0;

    EXPECT_EQ(1.0, cd[1]);
    EXPECT_EQ(2.0, cd[5]);
}

//---------------------------------------------------------------------------//

TEST_F(State_Test, 1PE_5Grp)
{
    if (nodes != 1)
        return;

    init(1, 1, 5);

    EXPECT_EQ(600, state->size());
    EXPECT_EQ(120, state->num_cells());
    EXPECT_EQ(5, state->num_groups());

    const State &const_state = *state;
    const_View_Field cd      = const_state.flux();
    for (int c = 0; c < 600; ++c)
    {
        EXPECT_EQ(0.0, cd[c]);
    }

    View_Field field = state->flux(1, 3);
    EXPECT_EQ(360, field.size());
    for (int g = 0; g < 3; ++g)
    {
        for (int c = 0; c < 120; ++c)
        {
            field[state->index(c, g)] = static_cast<double>(g+1);
        }
    }

    for (int c = 0; c < 120; ++c)
    {
        EXPECT_EQ(0.0, cd[state->index(c, 0)]);
        EXPECT_EQ(0.0, cd[state->index(c, 4)]);

        for (int g = 1; g <= 3; ++g)
        {
            EXPECT_EQ(static_cast<double>(g), cd[state->index(c, g)]);
        }
    }

    for (int g = 0; g < 5; ++g)
    {
        const_View_Field s = const_state.flux(g, g);
        EXPECT_EQ(120, s.size());

        for (int c = 0; c < 120; ++c)
        {
            EXPECT_EQ(c + 120 * g, state->index(c, g));
        }
    } 
}

//---------------------------------------------------------------------------//

TEST_F(State_Test, 4PE_5Grp)
{
    if (nodes != 4)
        return;

    init(2, 2, 5);

    EXPECT_EQ(5, state->num_groups());

    if (node == 0)
    {
        EXPECT_EQ(180, state->size());
        EXPECT_EQ(36, state->num_cells());
    }
    else if (node == 1)
    {
        EXPECT_EQ(120, state->size());
        EXPECT_EQ(24, state->num_cells());
    }
    else if (node == 2)
    {
        EXPECT_EQ(180, state->size());
        EXPECT_EQ(36, state->num_cells());
    }
    else if (node == 3)
    {
        EXPECT_EQ(120, state->size());
        EXPECT_EQ(24, state->num_cells());
    }
}

//---------------------------------------------------------------------------//
//                 end of tstState.cc
//---------------------------------------------------------------------------//
