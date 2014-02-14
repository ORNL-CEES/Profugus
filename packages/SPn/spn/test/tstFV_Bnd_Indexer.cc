//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstFV_Bnd_Indexer.cc
 * \author Thomas M. Evans
 * \date   Sat Nov 24 14:30:59 2012
 * \brief  FV_Bnd_Indexer unit test.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "utils/Definitions.hh"
#include "mesh/Partitioner.hh"
#include "../FV_Bnd_Indexer.hh"

using profugus::FV_Bnd_Indexer;

using namespace std;

using def::I; using def::J; using def::K;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class FV_Bnd_IndexerTest : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture

    typedef profugus::Partitioner          Partitioner;
    typedef Partitioner::ParameterList     ParameterList;
    typedef Partitioner::RCP_ParameterList RCP_ParameterList;
    typedef Partitioner::RCP_Mesh          RCP_Mesh;
    typedef Partitioner::RCP_Indexer       RCP_Indexer;
    typedef Partitioner::RCP_Global_Data   RCP_Global_Data;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        // build 4x4x4 mesh
        db = Teuchos::rcp(new ParameterList("test"));
        db->set("delta_x", 1.0);
        db->set("delta_y", 1.0);
        db->set("delta_z", 1.0);

        db->set("num_cells_i", 4);
        db->set("num_cells_j", 3);
        db->set("num_cells_k", 2);

        db->set("num_groups", 3);

        if (nodes == 2)
        {
            db->set("num_blocks_i", 2);
        }
        if (nodes == 4)
        {
            db->set("num_blocks_i", 2);
            db->set("num_blocks_j", 2);
        }

        Partitioner p(db);
        p.build();

        mesh    = p.get_mesh();
        indexer = p.get_indexer();
        data    = p.get_global_data();
    }

  protected:

    RCP_Mesh        mesh;
    RCP_Indexer     indexer;
    RCP_Global_Data data;

    RCP_ParameterList db;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(FV_Bnd_IndexerTest, 1_Block)
{
    if (nodes != 1) return;

    // assume high-end vacuum b.c's
    int bc_local[6] = {0, 6, 0, 8, 0, 12};

    EXPECT_EQ(0, indexer->offset(I));
    EXPECT_EQ(0, indexer->offset(J));

    // make indexer on each face
    FV_Bnd_Indexer hix(1, 3, 2, bc_local, 3, 2, bc_local, 0, 0);
    FV_Bnd_Indexer hiy(3, 4, 2, bc_local, 4, 2, bc_local, 0, 0);
    FV_Bnd_Indexer hiz(5, 4, 3, bc_local, 4, 3, bc_local, 0, 0);

    EXPECT_EQ(6, hix.num_local());
    EXPECT_EQ(6, hix.num_global());

    EXPECT_EQ(8, hiy.num_local());
    EXPECT_EQ(8, hiy.num_global());

    EXPECT_EQ(12, hiz.num_local());
    EXPECT_EQ(12, hiz.num_global());

    int ctr = 0;

    for (int o = 0; o < 2; ++o)
    {
        for (int a = 0; a < 3; ++a)
        {
            EXPECT_EQ(hix.local(a, o), hix.global(a, o));
            EXPECT_EQ(hix.global(a, o), hix.l2g(a, o));
            EXPECT_EQ(ctr, hix.local(a, o));
            ++ctr;
        }
    }

    for (int o = 0; o < 2; ++o)
    {
        for (int a = 0; a < 4; ++a)
        {
            EXPECT_EQ(hiy.local(a, o), hiy.global(a, o));
            EXPECT_EQ(hiy.global(a, o), hiy.l2g(a, o));
            EXPECT_EQ(ctr, hiy.local(a, o));
            ++ctr;
        }
    }

    for (int o = 0; o < 3; ++o)
    {
        for (int a = 0; a < 4; ++a)
        {
            EXPECT_EQ(hiz.local(a, o), hiz.global(a, o));
            EXPECT_EQ(hiz.global(a, o), hiz.l2g(a, o));
            EXPECT_EQ(ctr, hiz.local(a, o));
            ++ctr;
        }
    }

    EXPECT_EQ(26, ctr);
}

//---------------------------------------------------------------------------//

TEST_F(FV_Bnd_IndexerTest, 2_Block)
{
    if (nodes != 2) return;

    // assume high-end vacuum b.c's
    int bc_global[6] = {0, 6, 0, 8, 0, 12};
    int bc_local[6]  = {0, 0, 0, 0, 0, 0};

    int ctr = 0;

    if (node == 0)
    {
        bc_local[3] = 4;
        bc_local[5] = 6;

        FV_Bnd_Indexer hiy(3, 2, 2, bc_local, 4, 2, bc_global, 0, 0);
        FV_Bnd_Indexer hiz(5, 2, 3, bc_local, 4, 3, bc_global, 0, 0);

        EXPECT_EQ(4, hiy.num_local());
        EXPECT_EQ(8, hiy.num_global());

        EXPECT_EQ(6, hiz.num_local());
        EXPECT_EQ(12, hiz.num_global());

        for (int o = 0; o < 2; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                EXPECT_EQ(ctr, hiy.local(a, o));
                ++ctr;

                EXPECT_EQ(hiy.l2g(a, o), hiy.global(a, o));
            }
        }
        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                EXPECT_EQ(ctr, hiz.local(a, o));
                ++ctr;

                EXPECT_EQ(hiz.l2g(a, o), hiz.global(a, o));
            }
        }
        EXPECT_EQ(10, ctr);

        EXPECT_EQ(6,  hiy.l2g(0, 0));
        EXPECT_EQ(7,  hiy.l2g(1, 0));
        EXPECT_EQ(10, hiy.l2g(0, 1));
        EXPECT_EQ(11, hiy.l2g(1, 1));

        EXPECT_EQ(14, hiz.l2g(0, 0));
        EXPECT_EQ(15, hiz.l2g(1, 0));
        EXPECT_EQ(18, hiz.l2g(0, 1));
        EXPECT_EQ(19, hiz.l2g(1, 1));
        EXPECT_EQ(22, hiz.l2g(0, 2));
        EXPECT_EQ(23, hiz.l2g(1, 2));
    }

    if (node == 1)
    {
        bc_local[1] = 6;
        bc_local[3] = 4;
        bc_local[5] = 6;

        FV_Bnd_Indexer hix(1, 3, 2, bc_local, 3, 2, bc_global, 0, 0);
        FV_Bnd_Indexer hiy(3, 2, 2, bc_local, 4, 2, bc_global, 2, 0);
        FV_Bnd_Indexer hiz(5, 2, 3, bc_local, 4, 3, bc_global, 2, 0);

        EXPECT_EQ(6, hix.num_local());
        EXPECT_EQ(6, hix.num_global());

        EXPECT_EQ(4, hiy.num_local());
        EXPECT_EQ(8, hiy.num_global());

        EXPECT_EQ(6, hiz.num_local());
        EXPECT_EQ(12, hiz.num_global());

        for (int o = 0; o < 2; ++o)
        {
            for (int a = 0; a < 3; ++a)
            {
                EXPECT_EQ(ctr, hix.local(a, o));
                ++ctr;

                EXPECT_EQ(hix.l2g(a, o), hix.global(a, o));
            }
        }
        for (int o = 0; o < 2; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                EXPECT_EQ(ctr, hiy.local(a, o));
                ++ctr;

                EXPECT_EQ(hiy.l2g(a, o), hiy.global(a + indexer->offset(I), o));
            }
        }
        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                EXPECT_EQ(ctr, hiz.local(a, o));
                ++ctr;

                EXPECT_EQ(hiz.l2g(a, o), hiz.global(a + indexer->offset(I), o));
            }
        }
        EXPECT_EQ(16, ctr);

        EXPECT_EQ(0,  hix.l2g(0, 0));
        EXPECT_EQ(1,  hix.l2g(1, 0));
        EXPECT_EQ(2,  hix.l2g(2, 0));
        EXPECT_EQ(3,  hix.l2g(0, 1));
        EXPECT_EQ(4,  hix.l2g(1, 1));
        EXPECT_EQ(5,  hix.l2g(2, 1));

        EXPECT_EQ(8,  hiy.l2g(0, 0));
        EXPECT_EQ(9,  hiy.l2g(1, 0));
        EXPECT_EQ(12, hiy.l2g(0, 1));
        EXPECT_EQ(13, hiy.l2g(1, 1));

        EXPECT_EQ(16, hiz.l2g(0, 0));
        EXPECT_EQ(17, hiz.l2g(1, 0));
        EXPECT_EQ(20, hiz.l2g(0, 1));
        EXPECT_EQ(21, hiz.l2g(1, 1));
        EXPECT_EQ(24, hiz.l2g(0, 2));
        EXPECT_EQ(25, hiz.l2g(1, 2));
    }
}

//---------------------------------------------------------------------------//
//                        end of tstFV_Bnd_Indexer.cc
//---------------------------------------------------------------------------//
