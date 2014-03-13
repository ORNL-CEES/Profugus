//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstPartitioner.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 12 09:55:23 2014
 * \brief  Partitioner unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>
#include <cmath>
#include <sstream>

#include "Teuchos_RCP.hpp"

#include "utils/Definitions.hh"
#include "../Partitioner.hh"

using namespace std;

using def::I;
using def::J;
using def::K;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class Partitioner_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::Partitioner             Partitioner;
    typedef Teuchos::RCP<Partitioner>         RCP_Partitioner;
    typedef Partitioner::RCP_Mesh             RCP_Mesh;
    typedef Partitioner::RCP_Indexer          RCP_Indexer;
    typedef Partitioner::RCP_Global_Data      RCP_Global_Data;
    typedef Partitioner::RCP_ParameterList    RCP_ParameterList;
    typedef Partitioner::ParameterList        ParameterList;
    typedef Partitioner::Mesh_t::Space_Vector Space_Vector;
    typedef Partitioner::Mesh_t               Mesh;
    typedef Partitioner::Indexer_t            LG_Indexer;
    typedef Partitioner::Global_Data_t        Global_Mesh_Data;
    typedef Partitioner::Array_Dbl            Array_Dbl;


  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        pl = Teuchos::rcp(new ParameterList("Part"));
    }

  protected:
    // >>> Data that get re-initialized between tests

    RCP_Partitioner p;

    RCP_ParameterList pl;
    RCP_Mesh          mesh;
    RCP_Indexer       indexer;
    RCP_Global_Data   gdata;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// TESTS
//
//   j
//       _____
//       |   |  global-mesh    : 5 x 8 x 3
//       |   |  each partition : 5 x 8 x 3
//   0   |   |  dx       = 0.1
//       |   |  dy       = 0.2
//       |___|  dz       = 0.1
//              z-blocks = 1
//   i     0
//
TEST_F(Partitioner_Test, 1PE)
{
    if (nodes != 1)
        return;

    // set data
    {
        pl->set("num_blocks_i", 1);
        pl->set("num_blocks_j", 1);
        pl->set("num_cells_i", 5);
        pl->set("num_cells_j", 8);
        pl->set("num_cells_k", 3);
        pl->set("delta_x", 0.1);
        pl->set("delta_y", 0.2);
        pl->set("delta_z", 0.1);
        pl->set("num_z_blocks", 1);
    }

    // make the simple partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh
    mesh = p->get_mesh();
    EXPECT_FALSE(mesh.is_null());

    // check the mesh
    EXPECT_EQ(120, mesh->num_cells());
    EXPECT_EQ(5, mesh->num_cells_dim(I));
    EXPECT_EQ(8, mesh->num_cells_dim(J));
    EXPECT_EQ(3, mesh->num_cells_dim(K));

    EXPECT_EQ(0, mesh->block(def::I));
    EXPECT_EQ(0, mesh->block(def::J));
    EXPECT_EQ(1, mesh->block(def::K));

    EXPECT_EQ(5, mesh->num_cells_block_dim(I));
    EXPECT_EQ(8, mesh->num_cells_block_dim(J));
    EXPECT_EQ(3, mesh->num_cells_block_dim(K));

    for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
            {
                EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
            }
        }
    }

    Space_Vector block = mesh->block_widths();

    EXPECT_TRUE(soft_equiv(block[I], 0.5));
    EXPECT_TRUE(soft_equiv(block[J], 1.6));
    EXPECT_TRUE(soft_equiv(block[K], 0.3));

    EXPECT_TRUE(soft_equiv(mesh->center(0, I), 0.05));
    EXPECT_TRUE(soft_equiv(mesh->center(1, I), 0.15));
    EXPECT_TRUE(soft_equiv(mesh->center(2, I), 0.25));
    EXPECT_TRUE(soft_equiv(mesh->center(3, I), 0.35));
    EXPECT_TRUE(soft_equiv(mesh->center(4, I), 0.45));

    EXPECT_TRUE(soft_equiv(mesh->center(0, J), 0.1));
    EXPECT_TRUE(soft_equiv(mesh->center(1, J), 0.3));
    EXPECT_TRUE(soft_equiv(mesh->center(2, J), 0.5));
    EXPECT_TRUE(soft_equiv(mesh->center(3, J), 0.7));
    EXPECT_TRUE(soft_equiv(mesh->center(4, J), 0.9));
    EXPECT_TRUE(soft_equiv(mesh->center(5, J), 1.1));
    EXPECT_TRUE(soft_equiv(mesh->center(6, J), 1.3));
    EXPECT_TRUE(soft_equiv(mesh->center(7, J), 1.5));

    EXPECT_TRUE(soft_equiv(mesh->center(0, K), 0.05));
    EXPECT_TRUE(soft_equiv(mesh->center(1, K), 0.15));
    EXPECT_TRUE(soft_equiv(mesh->center(2, K), 0.25));
}

//---------------------------------------------------------------------------//
//
//   j
//       _________
//       |   |   |  global-mesh    : 10 x 8 x 3
//       |   |   |  each partition :  5 x 8 x 3
//   0   |   |   |  dx       = 0.1
//       |   |   |  dy       = 0.2
//       |___|___|  dz       = 0.1
//                  z-blocks = 1
//   i     0   1
//
TEST_F(Partitioner_Test, 2PE)
{
    if (nodes != 2)
        return;

    // set data
    {
        pl->set("num_blocks_i", 2);
        pl->set("num_blocks_j", 1);
        pl->set("num_cells_i", 10);
        pl->set("num_cells_j", 8);
        pl->set("num_cells_k", 3);
        pl->set("delta_x", 0.1);
        pl->set("delta_y", 0.2);
        pl->set("delta_z", 0.1);
        pl->set("num_z_blocks", 1);
    }

    // make the simple partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh
    mesh = p->get_mesh();
    EXPECT_FALSE(mesh.is_null());

    // check the mesh
    if (node == 0)
    {
        EXPECT_EQ(120, mesh->num_cells());
        EXPECT_EQ(5, mesh->num_cells_dim(I));
        EXPECT_EQ(8, mesh->num_cells_dim(J));
        EXPECT_EQ(3, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(5, mesh->num_cells_block_dim(I));
        EXPECT_EQ(8, mesh->num_cells_block_dim(J));
        EXPECT_EQ(3, mesh->num_cells_block_dim(K));

        for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
                {
                    EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                    EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                    EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
                }
            }
        }

        Space_Vector block = mesh->block_widths();

        EXPECT_TRUE(soft_equiv(block[I], 0.5));
        EXPECT_TRUE(soft_equiv(block[J], 1.6));
        EXPECT_TRUE(soft_equiv(block[K], 0.3));

        EXPECT_TRUE(soft_equiv(mesh->center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(mesh->center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(mesh->center(2, I), 0.25));
        EXPECT_TRUE(soft_equiv(mesh->center(3, I), 0.35));
        EXPECT_TRUE(soft_equiv(mesh->center(4, I), 0.45));

        EXPECT_TRUE(soft_equiv(mesh->center(0, J), 0.1));
        EXPECT_TRUE(soft_equiv(mesh->center(1, J), 0.3));
        EXPECT_TRUE(soft_equiv(mesh->center(2, J), 0.5));
        EXPECT_TRUE(soft_equiv(mesh->center(3, J), 0.7));
        EXPECT_TRUE(soft_equiv(mesh->center(4, J), 0.9));
        EXPECT_TRUE(soft_equiv(mesh->center(5, J), 1.1));
        EXPECT_TRUE(soft_equiv(mesh->center(6, J), 1.3));
        EXPECT_TRUE(soft_equiv(mesh->center(7, J), 1.5));

        EXPECT_TRUE(soft_equiv(mesh->center(0, K), 0.05));
        EXPECT_TRUE(soft_equiv(mesh->center(1, K), 0.15));
        EXPECT_TRUE(soft_equiv(mesh->center(2, K), 0.25));
    }
    else
    {
        EXPECT_EQ(120, mesh->num_cells());
        EXPECT_EQ(5, mesh->num_cells_dim(I));
        EXPECT_EQ(8, mesh->num_cells_dim(J));
        EXPECT_EQ(3, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(5, mesh->num_cells_block_dim(I));
        EXPECT_EQ(8, mesh->num_cells_block_dim(J));
        EXPECT_EQ(3, mesh->num_cells_block_dim(K));

        for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
                {
                    EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                    EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                    EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
                }
            }
        }

        Space_Vector block = mesh->block_widths();

        EXPECT_TRUE(soft_equiv(block[I], 0.5));
        EXPECT_TRUE(soft_equiv(block[J], 1.6));
        EXPECT_TRUE(soft_equiv(block[K], 0.3));

        EXPECT_TRUE(soft_equiv(mesh->center(0, I), 0.55));
        EXPECT_TRUE(soft_equiv(mesh->center(1, I), 0.65));
        EXPECT_TRUE(soft_equiv(mesh->center(2, I), 0.75));
        EXPECT_TRUE(soft_equiv(mesh->center(3, I), 0.85));
        EXPECT_TRUE(soft_equiv(mesh->center(4, I), 0.95));

        EXPECT_TRUE(soft_equiv(mesh->center(0, J), 0.1));
        EXPECT_TRUE(soft_equiv(mesh->center(1, J), 0.3));
        EXPECT_TRUE(soft_equiv(mesh->center(2, J), 0.5));
        EXPECT_TRUE(soft_equiv(mesh->center(3, J), 0.7));
        EXPECT_TRUE(soft_equiv(mesh->center(4, J), 0.9));
        EXPECT_TRUE(soft_equiv(mesh->center(5, J), 1.1));
        EXPECT_TRUE(soft_equiv(mesh->center(6, J), 1.3));
        EXPECT_TRUE(soft_equiv(mesh->center(7, J), 1.5));

        EXPECT_TRUE(soft_equiv(mesh->center(0, K), 0.05));
        EXPECT_TRUE(soft_equiv(mesh->center(1, K), 0.15));
        EXPECT_TRUE(soft_equiv(mesh->center(2, K), 0.25));
    }
}

//---------------------------------------------------------------------------//

TEST_F(Partitioner_Test, 4PE)
{
    if (nodes != 4)
        return;

    // set data
    {
        pl->set("num_blocks_i", 2);
        pl->set("num_blocks_j", 2);
        pl->set("num_cells_i", 4);
        pl->set("num_cells_j", 3);
        pl->set("num_cells_k", 2);
        pl->set("delta_x", 0.1);
        pl->set("delta_y", 0.2);
        pl->set("delta_z", 0.1);
        pl->set("num_z_blocks", 1);
    }

    // make the partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh
    mesh = p->get_mesh();
    EXPECT_FALSE(mesh.is_null());

    const Mesh &m = *mesh;

    // checks
    if (node == 0)
    {
        EXPECT_EQ(8, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
                {
                    EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                    EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                    EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
                }
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.4));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.2));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 0.3));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.15));
    }
    else if (node == 1)
    {
        EXPECT_EQ(8, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
                {
                    EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                    EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                    EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
                }
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.4));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.2));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.35));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 0.3));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.15));
    }
    else if (node == 2)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
                {
                    EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                    EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                    EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
                }
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.2));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.5));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.15));
    }
    else if (node == 3)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        for (int k = 0; k < mesh->num_cells_block_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
                {
                    EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                    EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
                    EXPECT_TRUE(soft_equiv(mesh->width(k, K), 0.1));
                }
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.2));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.35));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.5));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.15));
    }
    else
    {
        ADD_FAILURE();
    }

    // test the indexer
    indexer = p->get_indexer();
    EXPECT_FALSE(indexer.is_null());

    const LG_Indexer &i = *indexer;

    if (node == 0)
    {
        EXPECT_EQ(0, i.l2g(0, 0, 0));
        EXPECT_EQ(1, i.l2g(1, 0, 0));
        EXPECT_EQ(4, i.l2g(0, 1, 0));
        EXPECT_EQ(5, i.l2g(1, 1, 0));
        EXPECT_EQ(12, i.l2g(0, 0, 1));
        EXPECT_EQ(13, i.l2g(1, 0, 1));
        EXPECT_EQ(16, i.l2g(0, 1, 1));
        EXPECT_EQ(17, i.l2g(1, 1, 1));
    }
    else if (node == 1)
    {
        EXPECT_EQ(2, i.l2g(0, 0, 0));
        EXPECT_EQ(3, i.l2g(1, 0, 0));
        EXPECT_EQ(6, i.l2g(0, 1, 0));
        EXPECT_EQ(7, i.l2g(1, 1, 0));
        EXPECT_EQ(14, i.l2g(0, 0, 1));
        EXPECT_EQ(15, i.l2g(1, 0, 1));
        EXPECT_EQ(18, i.l2g(0, 1, 1));
        EXPECT_EQ(19, i.l2g(1, 1, 1));
    }
    else if (node == 2)
    {
        EXPECT_EQ(8, i.l2g(0, 0, 0));
        EXPECT_EQ(9, i.l2g(1, 0, 0));
        EXPECT_EQ(20, i.l2g(0, 0, 1));
        EXPECT_EQ(21, i.l2g(1, 0, 1));
    }
    else if (node == 3)
    {
        EXPECT_EQ(10, i.l2g(0, 0, 0));
        EXPECT_EQ(11, i.l2g(1, 0, 0));
        EXPECT_EQ(22, i.l2g(0, 0, 1));
        EXPECT_EQ(23, i.l2g(1, 0, 1));
    }
}

//---------------------------------------------------------------------------//

TEST_F(Partitioner_Test, 4PE_Nonuniform)
{
    if (nodes != 4)
        return;

    // set data
    {
        pl->set("num_blocks_i", 2);
        pl->set("num_blocks_j", 2);

        Array_Dbl x(5, 0.0);
        Array_Dbl y(4, 0.0);
        Array_Dbl z(3, 0.0);

        x[0] = -1.0;
        x[1] = 0.0;
        x[2] = 2.0;
        x[3] = 2.5;
        x[4] = 3.5;

        y[0] = 0.0;
        y[1] = 0.5;
        y[2] = 1.5;
        y[3] = 3.0;

        z[0] = -2.0;
        z[1] = 0.0;
        z[2] = 1.0;

        pl->set("x_edges", x);
        pl->set("y_edges", y);
        pl->set("z_edges", z);

        pl->set("num_z_blocks", 1);
    }

    // make the partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh
    mesh = p->get_mesh();
    EXPECT_FALSE(mesh.is_null());

    const Mesh &m = *mesh;

    // checks
    if (node == 0)
    {
        EXPECT_EQ(8, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 3.0));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 3.0));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), -0.5));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 1.0));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.center(0, K), -1.0));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.5));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 0.5));
        EXPECT_TRUE(soft_equiv(m.width(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 1.0));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), -1.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 0.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), -2.0));
    }
    else if (node == 1)
    {
        EXPECT_EQ(8, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 3.0));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 2.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 3.0));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.center(0, K), -1.0));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.5));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 0.5));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 0.5));
        EXPECT_TRUE(soft_equiv(m.width(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 1.0));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), 2.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 0.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), -2.0));
    }
    else if (node == 2)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 3.0));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 3.0));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), -0.5));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 1.0));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 2.25));
        EXPECT_TRUE(soft_equiv(m.center(0, K), -1.0));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.5));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 1.5));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 1.0));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), -1.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 1.5));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), -2.0));
    }
    else if (node == 3)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 3.0));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 2.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 3.00));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 2.25));
        EXPECT_TRUE(soft_equiv(m.center(0, K), -1.0));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.50));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 0.5));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 1.5));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 1.0));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), 2.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 1.5));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), -2.0));
    }
    else
    {
        ADD_FAILURE();
    }

    // test the indexer
    indexer = p->get_indexer();
    EXPECT_FALSE(indexer.is_null());

    const LG_Indexer &i = *indexer;

    if (node == 0)
    {
        EXPECT_EQ(0, i.l2g(0, 0, 0));
        EXPECT_EQ(1, i.l2g(1, 0, 0));
        EXPECT_EQ(4, i.l2g(0, 1, 0));
        EXPECT_EQ(5, i.l2g(1, 1, 0));
        EXPECT_EQ(12, i.l2g(0, 0, 1));
        EXPECT_EQ(13, i.l2g(1, 0, 1));
        EXPECT_EQ(16, i.l2g(0, 1, 1));
        EXPECT_EQ(17, i.l2g(1, 1, 1));
    }
    else if (node == 1)
    {
        EXPECT_EQ(2, i.l2g(0, 0, 0));
        EXPECT_EQ(3, i.l2g(1, 0, 0));
        EXPECT_EQ(6, i.l2g(0, 1, 0));
        EXPECT_EQ(7, i.l2g(1, 1, 0));
        EXPECT_EQ(14, i.l2g(0, 0, 1));
        EXPECT_EQ(15, i.l2g(1, 0, 1));
        EXPECT_EQ(18, i.l2g(0, 1, 1));
        EXPECT_EQ(19, i.l2g(1, 1, 1));
    }
    else if (node == 2)
    {
        EXPECT_EQ(8, i.l2g(0, 0, 0));
        EXPECT_EQ(9, i.l2g(1, 0, 0));
        EXPECT_EQ(20, i.l2g(0, 0, 1));
        EXPECT_EQ(21, i.l2g(1, 0, 1));
    }
    else if (node == 3)
    {
        EXPECT_EQ(10, i.l2g(0, 0, 0));
        EXPECT_EQ(11, i.l2g(1, 0, 0));
        EXPECT_EQ(22, i.l2g(0, 0, 1));
        EXPECT_EQ(23, i.l2g(1, 0, 1));
    }
}

//---------------------------------------------------------------------------//

TEST_F(Partitioner_Test, 4PE_Mixed)
{
    if (nodes != 4)
        return;

    // set data
    {
        pl->set("num_blocks_i", 2);
        pl->set("num_blocks_j", 2);

        Array_Dbl y(4, 0.0);

        y[0] = 0.0;
        y[1] = 0.5;
        y[2] = 1.5;
        y[3] = 3.0;

        pl->set("y_edges", y);

        pl->set("num_cells_i", 4);
        pl->set("num_cells_k", 2);

        pl->set("delta_x", 0.1);
        pl->set("delta_z", 0.2);

        pl->set("num_z_blocks", 1);
    }

    // make the partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh
    mesh = p->get_mesh();
    EXPECT_FALSE(mesh.is_null());

    const Mesh &m = *mesh;

    // checks
    if (node == 0)
    {
        EXPECT_EQ(8, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.4));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.3));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 0.5));
        EXPECT_TRUE(soft_equiv(m.width(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 0.2));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 0.2));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), 0.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 0.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), 0.0));
    }
    else if (node == 1)
    {
        EXPECT_EQ(8, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.4));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.35));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.3));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 0.5));
        EXPECT_TRUE(soft_equiv(m.width(1, J), 1.0));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 0.2));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 0.2));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), 0.2));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 0.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), 0.0));
    }
    else if (node == 2)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.4));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 2.25));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.3));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 1.5));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 0.2));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 0.2));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), 0.0));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 1.5));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), 0.0));
    }
    else if (node == 3)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));
        EXPECT_EQ(2, mesh->num_cells_dim(K));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));
        EXPECT_EQ(1, mesh->block(def::K));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));
        EXPECT_EQ(2, mesh->num_cells_block_dim(K));

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 1.5));
        EXPECT_TRUE(soft_equiv(bsize[K], 0.4));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.35));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 2.25));
        EXPECT_TRUE(soft_equiv(m.center(0, K), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, K), 0.3));

        // widths
        EXPECT_TRUE(soft_equiv(m.width(0, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(1, I), 0.1));
        EXPECT_TRUE(soft_equiv(m.width(0, J), 1.5));
        EXPECT_TRUE(soft_equiv(m.width(0, K), 0.2));
        EXPECT_TRUE(soft_equiv(m.width(1, K), 0.2));

        // corner
        EXPECT_TRUE(soft_equiv(m.low_corner(I), 0.2));
        EXPECT_TRUE(soft_equiv(m.low_corner(J), 1.5));
        EXPECT_TRUE(soft_equiv(m.low_corner(K), 0.0));
    }
    else
    {
        ADD_FAILURE();
    }

    // test the indexer
    indexer = p->get_indexer();
    EXPECT_FALSE(indexer.is_null());

    const LG_Indexer &i = *indexer;

    if (node == 0)
    {
        EXPECT_EQ(0, i.l2g(0, 0, 0));
        EXPECT_EQ(1, i.l2g(1, 0, 0));
        EXPECT_EQ(4, i.l2g(0, 1, 0));
        EXPECT_EQ(5, i.l2g(1, 1, 0));
        EXPECT_EQ(12, i.l2g(0, 0, 1));
        EXPECT_EQ(13, i.l2g(1, 0, 1));
        EXPECT_EQ(16, i.l2g(0, 1, 1));
        EXPECT_EQ(17, i.l2g(1, 1, 1));
    }
    else if (node == 1)
    {
        EXPECT_EQ(2, i.l2g(0, 0, 0));
        EXPECT_EQ(3, i.l2g(1, 0, 0));
        EXPECT_EQ(6, i.l2g(0, 1, 0));
        EXPECT_EQ(7, i.l2g(1, 1, 0));
        EXPECT_EQ(14, i.l2g(0, 0, 1));
        EXPECT_EQ(15, i.l2g(1, 0, 1));
        EXPECT_EQ(18, i.l2g(0, 1, 1));
        EXPECT_EQ(19, i.l2g(1, 1, 1));
    }
    else if (node == 2)
    {
        EXPECT_EQ(8, i.l2g(0, 0, 0));
        EXPECT_EQ(9, i.l2g(1, 0, 0));
        EXPECT_EQ(20, i.l2g(0, 0, 1));
        EXPECT_EQ(21, i.l2g(1, 0, 1));
    }
    else if (node == 3)
    {
        EXPECT_EQ(10, i.l2g(0, 0, 0));
        EXPECT_EQ(11, i.l2g(1, 0, 0));
        EXPECT_EQ(22, i.l2g(0, 0, 1));
        EXPECT_EQ(23, i.l2g(1, 0, 1));
    }

    // test the global mesh data
    gdata = p->get_global_data();

    EXPECT_EQ(24, gdata->num_cells());
    EXPECT_EQ(60, gdata->num_vertices());
    EXPECT_EQ(4, gdata->num_cells(I));
    EXPECT_EQ(3, gdata->num_cells(J));
    EXPECT_EQ(2, gdata->num_cells(K));
    EXPECT_EQ(5, gdata->num_vertices(I));
    EXPECT_EQ(4, gdata->num_vertices(J));
    EXPECT_EQ(3, gdata->num_vertices(K));
    EXPECT_EQ(0.0, gdata->low_edge(I));
    EXPECT_EQ(0.0, gdata->low_edge(J));
    EXPECT_EQ(0.0, gdata->low_edge(K));
    EXPECT_EQ(0.4, gdata->high_edge(I));
    EXPECT_EQ(3.0, gdata->high_edge(J));
    EXPECT_EQ(0.4, gdata->high_edge(K));
    EXPECT_TRUE(soft_equiv(gdata->volume(), 0.48));

    vector<double> x = gdata->edges(I);
    vector<double> y = gdata->edges(J);
    vector<double> z = gdata->edges(K);
    vector<double> rx(5, 0.0), ry(4, 0.0), rz(3, 0.0);
    {
        rx[1] = 0.1; rx[2] = 0.2; rx[3] = 0.3, rx[4] = 0.4;
        ry[1] = 0.5; ry[2] = 1.5; ry[3] = 3.0;
        rz[1] = 0.2; rz[2] = 0.4;
    }
    EXPECT_TRUE(soft_equiv(x.begin(), x.end(), rx.begin(), rx.end()));
    EXPECT_TRUE(soft_equiv(y.begin(), y.end(), ry.begin(), ry.end()));
    EXPECT_TRUE(soft_equiv(z.begin(), z.end(), rz.begin(), rz.end()));
}

//---------------------------------------------------------------------------//

TEST_F(Partitioner_Test, 4Domain)
{
    if (nodes != 4)
        return;

    // set data
    {
        pl->set("num_blocks_i", 2);
        pl->set("num_blocks_j", 1);

        Array_Dbl y(4, 0.0);

        y[0] = 0.0;
        y[1] = 0.5;
        y[2] = 1.5;
        y[3] = 3.0;

        pl->set("y_edges", y);

        pl->set("num_cells_i", 3);
        pl->set("num_cells_k", 2);

        pl->set("delta_x", 0.1);
        pl->set("delta_z", 0.2);

        pl->set("num_z_blocks", 2);

        pl->set("num_sets", 2);
    }

    // make the partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh, indexer, and global mesh
    mesh    = p->get_mesh();
    indexer = p->get_indexer();
    gdata   = p->get_global_data();

    EXPECT_FALSE(mesh.is_null());
    EXPECT_FALSE(indexer.is_null());
    EXPECT_FALSE(gdata.is_null());

    const Mesh             &m = *mesh;
    const LG_Indexer       &i = *indexer;
    const Global_Mesh_Data &g = *gdata;

    // check domains
    EXPECT_EQ(2, i.num_sets());
    EXPECT_EQ(2, i.num_blocks());
    EXPECT_EQ(3, i.num_global(I));
    EXPECT_EQ(3, i.num_global(J));
    EXPECT_EQ(2, i.num_blocks(I));
    EXPECT_EQ(1, i.num_blocks(J));
    EXPECT_EQ(2, i.num_cells_per_block(I)[0]);
    EXPECT_EQ(1, i.num_cells_per_block(I)[1]);
    EXPECT_EQ(3, i.num_cells_per_block(J)[0]);
    EXPECT_EQ(node, i.current_domain());

    EXPECT_EQ(18, g.num_cells());
    EXPECT_EQ(3, g.num_cells(I));
    EXPECT_EQ(3, g.num_cells(J));
    EXPECT_EQ(2, g.num_cells(K));

    if (node == 0)
    {
        EXPECT_EQ(0, i.set());
        EXPECT_EQ(0, i.block());

        EXPECT_EQ(0, i.l2g(0,0,0));
        EXPECT_EQ(1, i.l2g(1,0,0));
        EXPECT_EQ(3, i.l2g(0,1,0));
        EXPECT_EQ(4, i.l2g(1,1,0));
        EXPECT_EQ(6, i.l2g(0,2,0));
        EXPECT_EQ(7, i.l2g(1,2,0));
        EXPECT_EQ(9, i.l2g(0,0,1));
        EXPECT_EQ(10, i.l2g(1,0,1));
        EXPECT_EQ(12, i.l2g(0,1,1));
        EXPECT_EQ(13, i.l2g(1,1,1));
        EXPECT_EQ(15, i.l2g(0,2,1));
        EXPECT_EQ(16, i.l2g(1,2,1));

        EXPECT_EQ(12, m.num_cells());
        EXPECT_EQ(0, m.block(def::I));
        EXPECT_EQ(0, m.block(def::J));
    }
    else if (node == 1)
    {
        EXPECT_EQ(0, i.set());
        EXPECT_EQ(1, i.block());

        EXPECT_EQ(2, i.l2g(0,0,0));
        EXPECT_EQ(5, i.l2g(0,1,0));
        EXPECT_EQ(8, i.l2g(0,2,0));
        EXPECT_EQ(11, i.l2g(0,0,1));
        EXPECT_EQ(14, i.l2g(0,1,1));
        EXPECT_EQ(17, i.l2g(0,2,1));

        EXPECT_EQ(6, m.num_cells());
        EXPECT_EQ(1, m.block(def::I));
        EXPECT_EQ(0, m.block(def::J));
    }
    else if (node == 2)
    {
        EXPECT_EQ(1, i.set());
        EXPECT_EQ(0, i.block());

        EXPECT_EQ(0, i.l2g(0,0,0));
        EXPECT_EQ(1, i.l2g(1,0,0));
        EXPECT_EQ(3, i.l2g(0,1,0));
        EXPECT_EQ(4, i.l2g(1,1,0));
        EXPECT_EQ(6, i.l2g(0,2,0));
        EXPECT_EQ(7, i.l2g(1,2,0));
        EXPECT_EQ(9, i.l2g(0,0,1));
        EXPECT_EQ(10, i.l2g(1,0,1));
        EXPECT_EQ(12, i.l2g(0,1,1));
        EXPECT_EQ(13, i.l2g(1,1,1));
        EXPECT_EQ(15, i.l2g(0,2,1));
        EXPECT_EQ(16, i.l2g(1,2,1));

        EXPECT_EQ(12, m.num_cells());
        EXPECT_EQ(0, m.block(def::I));
        EXPECT_EQ(0, m.block(def::J));
    }
    else if (node == 3)
    {
        EXPECT_EQ(1, i.set());
        EXPECT_EQ(1, i.block());

        EXPECT_EQ(2, i.l2g(0,0,0));
        EXPECT_EQ(5, i.l2g(0,1,0));
        EXPECT_EQ(8, i.l2g(0,2,0));
        EXPECT_EQ(11, i.l2g(0,0,1));
        EXPECT_EQ(14, i.l2g(0,1,1));
        EXPECT_EQ(17, i.l2g(0,2,1));

        EXPECT_EQ(6, m.num_cells());
        EXPECT_EQ(1, m.block(def::I));
        EXPECT_EQ(0, m.block(def::J));
    }
}

//---------------------------------------------------------------------------//

TEST_F(Partitioner_Test, 4PE_2D)
{
    if (nodes != 4)
        return;

    // set data
    {
        pl->set("num_blocks_i", 2);
        pl->set("num_blocks_j", 2);
        pl->set("num_cells_i", 4);
        pl->set("num_cells_j", 3);
        pl->set("delta_x", 0.1);
        pl->set("delta_y", 0.2);
    }

    // make the partitioner
    {
        // simple partitioner specialization
        p = Teuchos::rcp(new Partitioner(pl));
        EXPECT_FALSE(p.is_null());
    }

    // build the mesh
    p->build();

    // get the mesh
    mesh = p->get_mesh();
    EXPECT_FALSE(mesh.is_null());

    const Mesh &m = *mesh;

    // checks
    if (node == 0)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));

        for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
            {
                EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.4));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 0.3));
    }
    else if (node == 1)
    {
        EXPECT_EQ(4, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(2, mesh->num_cells_dim(J));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(0, mesh->block(def::J));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(2, mesh->num_cells_block_dim(J));

        for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
            {
                EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.4));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.35));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.1));
        EXPECT_TRUE(soft_equiv(m.center(1, J), 0.3));
    }
    else if (node == 2)
    {
        EXPECT_EQ(2, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));

        EXPECT_EQ(0, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));

        for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
            {
                EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.2));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.05));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.15));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.5));
    }
    else if (node == 3)
    {
        EXPECT_EQ(2, mesh->num_cells());
        EXPECT_EQ(2, mesh->num_cells_dim(I));
        EXPECT_EQ(1, mesh->num_cells_dim(J));

        EXPECT_EQ(1, mesh->block(def::I));
        EXPECT_EQ(1, mesh->block(def::J));

        EXPECT_EQ(2, mesh->num_cells_block_dim(I));
        EXPECT_EQ(1, mesh->num_cells_block_dim(J));

        for (int j = 0; j < mesh->num_cells_block_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_block_dim(I); ++i)
            {
                EXPECT_TRUE(soft_equiv(mesh->width(i, I), 0.1));
                EXPECT_TRUE(soft_equiv(mesh->width(j, J), 0.2));
            }
        }

        // check the block size
        Space_Vector bsize = m.block_widths();

        EXPECT_TRUE(soft_equiv(bsize[I], 0.2));
        EXPECT_TRUE(soft_equiv(bsize[J], 0.2));

        // centers
        EXPECT_TRUE(soft_equiv(m.center(0, I), 0.25));
        EXPECT_TRUE(soft_equiv(m.center(1, I), 0.35));
        EXPECT_TRUE(soft_equiv(m.center(0, J), 0.5));
    }
    else
    {
        ADD_FAILURE();
    }

    // test the indexer
    indexer = p->get_indexer();
    EXPECT_FALSE(indexer.is_null());

    const LG_Indexer &i = *indexer;

    if (node == 0)
    {
        EXPECT_EQ(0, i.l2g(0, 0, 0));
        EXPECT_EQ(1, i.l2g(1, 0, 0));
        EXPECT_EQ(4, i.l2g(0, 1, 0));
        EXPECT_EQ(5, i.l2g(1, 1, 0));
    }
    else if (node == 1)
    {
        EXPECT_EQ(2, i.l2g(0, 0, 0));
        EXPECT_EQ(3, i.l2g(1, 0, 0));
        EXPECT_EQ(6, i.l2g(0, 1, 0));
        EXPECT_EQ(7, i.l2g(1, 1, 0));
    }
    else if (node == 2)
    {
        EXPECT_EQ(8, i.l2g(0, 0, 0));
        EXPECT_EQ(9, i.l2g(1, 0, 0));
    }
    else if (node == 3)
    {
        EXPECT_EQ(10, i.l2g(0, 0, 0));
        EXPECT_EQ(11, i.l2g(1, 0, 0));
    }
}

//---------------------------------------------------------------------------//
//                 end of tstPartitioner.cc
//---------------------------------------------------------------------------//
