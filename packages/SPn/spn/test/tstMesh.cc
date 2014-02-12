//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstMesh.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 12 00:11:40 2014
 * \brief  Mesh unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Teuchos_RCP.hpp"

#include "utils/Definitions.hh"
#include "../Mesh.hh"

using def::I;
using def::J;
using def::K;

//---------------------------------------------------------------------------//
// FIXTURE
//---------------------------------------------------------------------------//

class Mesh_Test : public testing::Test
{
  protected:
    typedef profugus::Mesh     Mesh;
    typedef Teuchos::RCP<Mesh> RCP_Mesh;

  protected:
    Mesh_Test()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

    virtual ~Mesh_Test()
    {
    }

  protected:

    int node, nodes;
};

//---------------------------------------------------------------------------//

class Mesh_3D : public Mesh_Test
{
};

//---------------------------------------------------------------------------//

class Mesh_2D : public Mesh_Test
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Mesh_3D, 1PE)
{
    if (nodes > 1)
        return;

    RCP_Mesh mesh;

    // make a 4x3x2 mesh
    {
        def::Vec_Dbl x(5, 0.0);
        def::Vec_Dbl y(4, 0.0);
        def::Vec_Dbl z(3, 0.0);

        x[1] = 1.0;
        x[2] = 2.0;
        x[3] = 3.0;
        x[4] = 4.0;

        y[1] = 2.0;
        y[2] = 3.0;
        y[3] = 3.5;

        z[1] = 0.5;
        z[2] = 1.0;

        mesh = Teuchos::rcp(new Mesh(x, y, z, 0, 0, 1));
    }

    const Mesh &m = *mesh;

    // check
    EXPECT_EQ(3, m.dimension());
    EXPECT_EQ(4, m.num_cells_dim(I));
    EXPECT_EQ(3, m.num_cells_dim(J));
    EXPECT_EQ(2, m.num_cells_dim(K));
    EXPECT_EQ(24, m.num_cells());
    EXPECT_EQ(60, m.num_vertices());
    EXPECT_EQ(0, m.block(I));
    EXPECT_EQ(0, m.block(J));
    EXPECT_EQ(1, m.block(K));
    EXPECT_EQ(4, m.num_cells_block_dim(I));
    EXPECT_EQ(3, m.num_cells_block_dim(J));
    EXPECT_EQ(2, m.num_cells_block_dim(K));

    unsigned long cell= 23;
    EXPECT_EQ(cell, m.convert(3,2,1));

    Mesh::Dim_Vector ijk = m.cardinal(cell);
    EXPECT_EQ(3, ijk[I]);
    EXPECT_EQ(2, ijk[J]);
    EXPECT_EQ(1, ijk[K]);
    EXPECT_TRUE(soft_equiv(m.volume(cell), 0.25));

    EXPECT_EQ("xyz_3d",  mesh->label());
}

//---------------------------------------------------------------------------//

TEST_F(Mesh_3D, 4PE)
{
    if (nodes != 4)
        return;

    RCP_Mesh mesh;

    // make a 2x2 mesh decomposition (in i,j) on 4x3x2 mesh
    if (node == 0)
    {
        def::Vec_Dbl x(3, 0.0);
        def::Vec_Dbl y(3, 0.0);
        def::Vec_Dbl z(3, 0.0);

        x[1] = 1.0;
        x[2] = 2.0;

        y[1] = 2.0;
        y[2] = 3.0;

        z[1] = 0.5;
        z[2] = 1.0;

        mesh = Teuchos::rcp(new Mesh(x, y, z, 0, 0, 2));
    }
    else if (node == 1)
    {
        def::Vec_Dbl x(3, 0.0);
        def::Vec_Dbl y(3, 0.0);
        def::Vec_Dbl z(3, 0.0);

        x[0] = 2.0;
        x[1] = 3.0;
        x[2] = 4.0;

        y[1] = 2.0;
        y[2] = 3.0;

        z[1] = 0.5;
        z[2] = 1.0;

        mesh = Teuchos::rcp(new Mesh(x, y, z, 0, 1, 2));
    }
    else if (node == 2)
    {
        def::Vec_Dbl x(3, 0.0);
        def::Vec_Dbl y(2, 0.0);
        def::Vec_Dbl z(3, 0.0);

        x[1] = 1.0;
        x[2] = 2.0;

        y[0] = 3.0;
        y[1] = 3.5;

        z[1] = 0.5;
        z[2] = 1.0;

        mesh = Teuchos::rcp(new Mesh(x, y, z, 1, 0, 2));
    }
    else if (node == 3)
   {
        def::Vec_Dbl x(3, 0.0);
        def::Vec_Dbl y(2, 0.0);
        def::Vec_Dbl z(3, 0.0);

        x[0] = 2.0;
        x[1] = 3.0;
        x[2] = 4.0;

        y[0] = 3.0;
        y[1] = 3.5;

        z[1] = 0.5;
        z[2] = 1.0;

        mesh = Teuchos::rcp(new Mesh(x, y, z, 1, 1, 2));
    }

    // check the mesh
    const Mesh &m = *mesh;

    EXPECT_EQ(2, m.block(K));

    if (node == 0)
    {
        EXPECT_EQ(2, m.num_cells_dim(I));
        EXPECT_EQ(2, m.num_cells_dim(J));
        EXPECT_EQ(2, m.num_cells_dim(K));
        EXPECT_EQ(8, m.num_cells());
        EXPECT_EQ(27, m.num_vertices());
        EXPECT_EQ(0, m.block(I));
        EXPECT_EQ(0, m.block(J));
        EXPECT_EQ(2, m.num_cells_block_dim(I));
        EXPECT_EQ(2, m.num_cells_block_dim(J));
        EXPECT_EQ(1, m.num_cells_block_dim(K));

        for (int i = 0; i < m.num_cells_dim(I); ++i)
            EXPECT_TRUE(soft_equiv(m.width(i, I), 1.0));

        EXPECT_TRUE(soft_equiv(m.width(0, J), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(1, J), 1.0));

        for (int k = 0; k < m.num_cells_dim(K); ++k)
            EXPECT_TRUE(soft_equiv(m.width(k, K), 0.5));

        // do some searches
        Mesh::Dim_Vector ijk(-1);

        EXPECT_TRUE(m.find_cell(Mesh::Space_Vector(1.5, 2.01, 0.01), ijk));
        EXPECT_TRUE(ijk[I] == 1);
        EXPECT_TRUE(ijk[J] == 1);
        EXPECT_TRUE(ijk[K] == 0);

        EXPECT_TRUE(!m.find_cell(Mesh::Space_Vector(2.01, 2.01, 0.01), ijk));

        EXPECT_TRUE(m.find_cell(Mesh::Space_Vector(1.5, 2.0, 0.6), ijk));
        EXPECT_TRUE(ijk[I] == 1);
        EXPECT_TRUE(ijk[J] == -1);
        EXPECT_TRUE(ijk[K] == 1);
    }
    else if (node == 1)
    {
        EXPECT_EQ(2, m.num_cells_dim(I));
        EXPECT_EQ(2, m.num_cells_dim(J));
        EXPECT_EQ(2, m.num_cells_dim(K));
        EXPECT_EQ(8, m.num_cells());
        EXPECT_EQ(0, m.block(I));
        EXPECT_EQ(1, m.block(J));
        EXPECT_EQ(2, m.num_cells_block_dim(I));
        EXPECT_EQ(2, m.num_cells_block_dim(J));
        EXPECT_EQ(1, m.num_cells_block_dim(K));

        for (int i = 0; i < m.num_cells_dim(I); ++i)
            EXPECT_TRUE(soft_equiv(m.width(i, I), 1.0));

        EXPECT_TRUE(soft_equiv(m.width(0, J), 2.0));
        EXPECT_TRUE(soft_equiv(m.width(1, J), 1.0));

        for (int k = 0; k < m.num_cells_dim(K); ++k)
            EXPECT_TRUE(soft_equiv(m.width(k, K), 0.5));

        EXPECT_TRUE(soft_equiv(m.center(1, J), 2.5));

        Mesh::Space_Vector bw = m.block_widths();
        EXPECT_TRUE(soft_equiv(bw[I], 2.0));
        EXPECT_TRUE(soft_equiv(bw[J], 3.0));
        EXPECT_TRUE(soft_equiv(bw[K], 1.0));

    }
    else if (node == 2)
    {
        EXPECT_EQ(2, m.num_cells_dim(I));
        EXPECT_EQ(1, m.num_cells_dim(J));
        EXPECT_EQ(2, m.num_cells_dim(K));
        EXPECT_EQ(4, m.num_cells());
        EXPECT_EQ(1, m.block(I));
        EXPECT_EQ(0, m.block(J));
        EXPECT_EQ(2, m.num_cells_block_dim(I));
        EXPECT_EQ(1, m.num_cells_block_dim(J));
        EXPECT_EQ(1, m.num_cells_block_dim(K));

        for (int i = 0; i < m.num_cells_dim(I); ++i)
            EXPECT_TRUE(soft_equiv(m.width(i, I), 1.0));

        for (int j = 0; j < m.num_cells_dim(J); ++j)
            EXPECT_TRUE(soft_equiv(m.width(j, J), 0.5));

        for (int k = 0; k < m.num_cells_dim(K); ++k)
            EXPECT_TRUE(soft_equiv(m.width(k, K), 0.5));

        // do some searches
        Mesh::Dim_Vector ijk(-1);

        EXPECT_TRUE(m.find_cell(Mesh::Space_Vector(0.5, 3.01, 1.0), ijk));
        EXPECT_TRUE(ijk[I] == 0);
        EXPECT_TRUE(ijk[J] == 0);
        EXPECT_TRUE(ijk[K] == -2);

        EXPECT_TRUE(!m.find_cell(Mesh::Space_Vector(2.01, 2.01, 0.01), ijk));

        EXPECT_TRUE(m.find_cell(Mesh::Space_Vector(1.5, 3.49, 0.6), ijk));
        EXPECT_TRUE(ijk[I] == 1);
        EXPECT_TRUE(ijk[J] == 0);
        EXPECT_TRUE(ijk[K] == 1);

        EXPECT_TRUE(m.find_cell(Mesh::Space_Vector(0.0, 3.0, 0.6), ijk));
        EXPECT_TRUE(ijk[I] == 0);
        EXPECT_TRUE(ijk[J] == 0);
        EXPECT_TRUE(ijk[K] == 1);
    }
   else if (node == 3)
    {
        EXPECT_EQ(2, m.num_cells_dim(I));
        EXPECT_EQ(1, m.num_cells_dim(J));
        EXPECT_EQ(2, m.num_cells_dim(K));
        EXPECT_EQ(4, m.num_cells());
        EXPECT_EQ(1, m.block(I));
        EXPECT_EQ(1, m.block(J));
        EXPECT_EQ(2, m.num_cells_block_dim(I));
        EXPECT_EQ(1, m.num_cells_block_dim(J));
        EXPECT_EQ(1, m.num_cells_block_dim(K));

        for (int i = 0; i < m.num_cells_dim(I); ++i)
            EXPECT_TRUE(soft_equiv(m.width(i, I), 1.0));
        for (int i = 0; i < m.num_cells_dim(I); ++i)
            EXPECT_TRUE(soft_equiv(m.inv_width(i, I), 1.0));

        for (int j = 0; j < m.num_cells_dim(J); ++j)
            EXPECT_TRUE(soft_equiv(m.width(j, J), 0.5));
        for (int j = 0; j < m.num_cells_dim(J); ++j)
            EXPECT_TRUE(soft_equiv(m.inv_width(j, J), 2.0));

        for (int k = 0; k < m.num_cells_dim(K); ++k)
            EXPECT_TRUE(soft_equiv(m.width(k, K), 0.5));
        for (int k = 0; k < m.num_cells_dim(K); ++k)
            EXPECT_TRUE(soft_equiv(m.inv_width(k, K), 2.0));
    }
}

//---------------------------------------------------------------------------//
//                 end of tstMesh.cc
//---------------------------------------------------------------------------//
