//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/test/tstSDM_Face_Field.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 30 23:04:27 2012
 * \brief  Unit-Test for SDM_Face_Field
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <algorithm>
#include <vector>

#include <SPn/config.h>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "mesh/Partitioner.hh"
#include "../SDM_Face_Field.hh"

using def::X;
using def::Y;
using def::Z;

using namespace std;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class SDM_Face_Field_Test : public testing::Test
{
  protected:
    typedef profugus::SDM_Face_Field       SDM_Face_Field;
    typedef profugus::Partitioner          Partitioner;
    typedef Partitioner::ParameterList     ParameterList;
    typedef Partitioner::RCP_ParameterList RCP_ParameterList;
    typedef Partitioner::RCP_Mesh          RCP_Mesh;
    typedef SDM_Face_Field::Serial_Matrix  Serial_Matrix;
    typedef Teuchos::RCP<SDM_Face_Field>   RCP_Face_Field;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        // build 4x4x4 mesh
        RCP_ParameterList db = Teuchos::rcp(new ParameterList("test"));

        db->set("delta_x", 1.0);
        db->set("delta_y", 1.0);
        db->set("delta_z", 1.0);

        db->set("num_cells_i", 6);
        db->set("num_cells_j", 4);
        db->set("num_cells_k", 2);

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

        mesh = p.get_mesh();

        x = Teuchos::rcp(new SDM_Face_Field(*mesh, X, 3));
        y = Teuchos::rcp(new SDM_Face_Field(*mesh, Y, 3));
        z = Teuchos::rcp(new SDM_Face_Field(*mesh, Z, 3));

        EXPECT_EQ(0, x->face());
        EXPECT_EQ(1, y->face());
        EXPECT_EQ(2, z->face());
    }

    void invert(Serial_Matrix &A)
    {
        EXPECT_EQ(A.numRows(), A.numCols());

        int G = A.numRows();

        // LAPACK work arrays.
        vector<double> work(G);
        vector<int>    ipiv(G);
        int info = 0;

        // LU decomposition
        lapack.GETRF(G, G, A.values(), A.stride(), &ipiv[0], &info);
        CHECK(info == 0);

        // inverse
        lapack.GETRI(G, A.values(), A.stride(), &ipiv[0], &work[0], G, &info);
        CHECK(info == 0);
    }

  protected:

    RCP_Mesh mesh;

    RCP_Face_Field x, y, z;

    int node, nodes;

    Teuchos::LAPACK<int, double> lapack;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SDM_Face_Field_Test, Basic_1pe)
{
    if (nodes != 1) return;

    EXPECT_EQ(8,  x->size());
    EXPECT_EQ(12, y->size());
    EXPECT_EQ(24, z->size());

    EXPECT_EQ(4, x->abscissa());
    EXPECT_EQ(2, x->ordinate());

    EXPECT_EQ(6, y->abscissa());
    EXPECT_EQ(2, y->ordinate());

    EXPECT_EQ(6, z->abscissa());
    EXPECT_EQ(4, z->ordinate());

    EXPECT_EQ(72,  x->data_size());
    EXPECT_EQ(108, y->data_size());
    EXPECT_EQ(216, z->data_size());
}

//---------------------------------------------------------------------------//

TEST_F(SDM_Face_Field_Test, Basic_2pe)
{
    if (nodes != 2) return;

    EXPECT_EQ(8,  x->size());
    EXPECT_EQ(6,  y->size());
    EXPECT_EQ(12, z->size());

    EXPECT_EQ(4, x->abscissa());
    EXPECT_EQ(2, x->ordinate());

    EXPECT_EQ(3, y->abscissa());
    EXPECT_EQ(2, y->ordinate());

    EXPECT_EQ(3, z->abscissa());
    EXPECT_EQ(4, z->ordinate());

    EXPECT_EQ(72,  x->data_size());
    EXPECT_EQ(54,  y->data_size());
    EXPECT_EQ(108, z->data_size());
}

//---------------------------------------------------------------------------//

TEST_F(SDM_Face_Field_Test, Basic_4pe)
{
    if (nodes != 4) return;

    EXPECT_EQ(4,  x->size());
    EXPECT_EQ(6,  y->size());
    EXPECT_EQ(6 , z->size());

    EXPECT_EQ(2, x->abscissa());
    EXPECT_EQ(2, x->ordinate());

    EXPECT_EQ(3, y->abscissa());
    EXPECT_EQ(2, y->ordinate());

    EXPECT_EQ(3, z->abscissa());
    EXPECT_EQ(2, z->ordinate());

    EXPECT_EQ(36, x->data_size());
    EXPECT_EQ(54, y->data_size());
    EXPECT_EQ(54, z->data_size());
}

//---------------------------------------------------------------------------//

TEST_F(SDM_Face_Field_Test, Data_Access)
{
    if (nodes != 1) return;

    EXPECT_EQ(4, x->abscissa());
    EXPECT_EQ(2, x->ordinate());

    // make some data
    {
        Serial_Matrix a(3, 3), b(3, 3), c(3, 3);
        a(0, 0) = 1.1; a(0, 1) = 0.0; a(0, 2) = 0.0;
        a(1, 0) = 1.4; a(1, 1) = 1.5; a(1, 2) = 1.6;
        a(2, 0) = 1.7; a(2, 1) = 1.8; a(2, 2) = 1.9;

        b(0, 0) = 2.1; b(0, 1) = 0.0; b(0, 2) = 0.0;
        b(1, 0) = 2.4; b(1, 1) = 2.5; b(1, 2) = 2.6;
        b(2, 0) = 2.7; b(2, 1) = 2.8; b(2, 2) = 2.9;

        c(0, 0) = 3.1; c(0, 1) = 0.0; c(0, 2) = 0.0;
        c(1, 0) = 3.4; c(1, 1) = 3.5; c(1, 2) = 3.6;
        c(2, 0) = 3.7; c(2, 1) = 3.8; c(2, 2) = 3.9;

        for (int k = 0; k < x->ordinate(); ++k)
        {
            for (int j = 0; j < x->abscissa() / 2; ++j)
            {
                x->insert(j, k, a);
            }
            for (int j = x->abscissa() / 2; j < x->abscissa(); ++j)
            {
                x->insert(j, k, b);
            }
        }

        // add c to the field
        x->insert(1, 1, c);
    }

    // check the data
    const double *data = x->data_pointer();
    for (SDM_Face_Field::const_pointer dp  = x->begin_data();
                                       dp != x->end_data(); ++dp)
    {
        EXPECT_EQ(*data, *dp);
        ++data;
    }

    data = x->data_pointer();

    // (0,0) starts at (0 + 0 * 4) * 9 = 0
    EXPECT_EQ(data[0], 1.1);
    EXPECT_EQ(data[1], 1.4);
    EXPECT_EQ(data[2], 1.7);
    EXPECT_EQ(data[3], 0.0);
    EXPECT_EQ(data[4], 1.5);
    EXPECT_EQ(data[5], 1.8);
    EXPECT_EQ(data[6], 0.0);
    EXPECT_EQ(data[7], 1.6);
    EXPECT_EQ(data[8], 1.9);

    // (2,0) starts at (2 + 0 * 4) * 9 = 18
    EXPECT_EQ(data[18], 2.1);
    EXPECT_EQ(data[19], 2.4);
    EXPECT_EQ(data[20], 2.7);
    EXPECT_EQ(data[21], 0.0);
    EXPECT_EQ(data[22], 2.5);
    EXPECT_EQ(data[23], 2.8);
    EXPECT_EQ(data[24], 0.0);
    EXPECT_EQ(data[25], 2.6);
    EXPECT_EQ(data[26], 2.9);

    // (0,1) starts at (0 + 1 * 4) * 9 = 36
    EXPECT_EQ(data[36], 1.1);
    EXPECT_EQ(data[37], 1.4);
    EXPECT_EQ(data[38], 1.7);
    EXPECT_EQ(data[39], 0.0);
    EXPECT_EQ(data[40], 1.5);
    EXPECT_EQ(data[41], 1.8);
    EXPECT_EQ(data[42], 0.0);
    EXPECT_EQ(data[43], 1.6);
    EXPECT_EQ(data[44], 1.9);

    // (1,1) starts at (1 + 1 * 4) * 9 = 45
    EXPECT_EQ(data[45], 3.1);
    EXPECT_EQ(data[46], 3.4);
    EXPECT_EQ(data[47], 3.7);
    EXPECT_EQ(data[48], 0.0);
    EXPECT_EQ(data[49], 3.5);
    EXPECT_EQ(data[50], 3.8);
    EXPECT_EQ(data[51], 0.0);
    EXPECT_EQ(data[52], 3.6);
    EXPECT_EQ(data[53], 3.9);

    // (2,1) starts at (2 + 1 * 4) * 9 = 54
    EXPECT_EQ(data[54], 2.1);
    EXPECT_EQ(data[55], 2.4);
    EXPECT_EQ(data[56], 2.7);
    EXPECT_EQ(data[57], 0.0);
    EXPECT_EQ(data[58], 2.5);
    EXPECT_EQ(data[59], 2.8);
    EXPECT_EQ(data[60], 0.0);
    EXPECT_EQ(data[61], 2.6);
    EXPECT_EQ(data[62], 2.9);

    // check matrix views
    {
        Serial_Matrix aa = x->view(0, 1);
        Serial_Matrix bb = x->view(2, 0);
        Serial_Matrix cc = x->view(1, 1);

        EXPECT_EQ(aa(0, 0), 1.1);
        EXPECT_EQ(aa(0, 1), 0.0);
        EXPECT_EQ(aa(0, 2), 0.0);
        EXPECT_EQ(aa(1, 0), 1.4);
        EXPECT_EQ(aa(1, 1), 1.5);
        EXPECT_EQ(aa(1, 2), 1.6);
        EXPECT_EQ(aa(2, 0), 1.7);
        EXPECT_EQ(aa(2, 1), 1.8);
        EXPECT_EQ(aa(2, 2), 1.9);

        EXPECT_EQ(bb(0, 0), 2.1);
        EXPECT_EQ(bb(0, 1), 0.0);
        EXPECT_EQ(bb(0, 2), 0.0);
        EXPECT_EQ(bb(1, 0), 2.4);
        EXPECT_EQ(bb(1, 1), 2.5);
        EXPECT_EQ(bb(1, 2), 2.6);
        EXPECT_EQ(bb(2, 0), 2.7);
        EXPECT_EQ(bb(2, 1), 2.8);
        EXPECT_EQ(bb(2, 2), 2.9);

        EXPECT_EQ(cc(0, 0), 3.1);
        EXPECT_EQ(cc(0, 1), 0.0);
        EXPECT_EQ(cc(0, 2), 0.0);
        EXPECT_EQ(cc(1, 0), 3.4);
        EXPECT_EQ(cc(1, 1), 3.5);
        EXPECT_EQ(cc(1, 2), 3.6);
        EXPECT_EQ(cc(2, 0), 3.7);
        EXPECT_EQ(cc(2, 1), 3.8);
        EXPECT_EQ(cc(2, 2), 3.9);

        // change the underlying data
        invert(bb);
    }

    // (0,0) starts at (0 + 0 * 4) * 9 = 0
    EXPECT_EQ(data[0], 1.1);
    EXPECT_EQ(data[1], 1.4);
    EXPECT_EQ(data[2], 1.7);
    EXPECT_EQ(data[3], 0.0);
    EXPECT_EQ(data[4], 1.5);
    EXPECT_EQ(data[5], 1.8);
    EXPECT_EQ(data[6], 0.0);
    EXPECT_EQ(data[7], 1.6);
    EXPECT_EQ(data[8], 1.9);

    double eps = 1.0e-10;

    // (2,0) starts at (2 + 0 * 4) * 9 = 18
    EXPECT_NEAR(data[18],  0.476190476190476, eps);
    EXPECT_NEAR(data[19], -0.952380952380973, eps);
    EXPECT_NEAR(data[20],  0.476190476190496, eps);
    EXPECT_NEAR(data[21],  0.0, eps);
    EXPECT_NEAR(data[22], -96.666666666667183, eps);
    EXPECT_NEAR(data[23],  93.333333333333826, eps);
    EXPECT_NEAR(data[24],  0.0, eps);
    EXPECT_NEAR(data[25],  86.666666666667140, eps);
    EXPECT_NEAR(data[26], -83.333333333333783, eps);

    // (0,1) starts at (0 + 1 * 4) * 9 = 36
    EXPECT_EQ(data[36], 1.1);
    EXPECT_EQ(data[37], 1.4);
    EXPECT_EQ(data[38], 1.7);
    EXPECT_EQ(data[39], 0.0);
    EXPECT_EQ(data[40], 1.5);
    EXPECT_EQ(data[41], 1.8);
    EXPECT_EQ(data[42], 0.0);
    EXPECT_EQ(data[43], 1.6);
    EXPECT_EQ(data[44], 1.9);

    // (1,1) starts at (1 + 1 * 4) * 9 = 45
    EXPECT_EQ(data[45], 3.1);
    EXPECT_EQ(data[46], 3.4);
    EXPECT_EQ(data[47], 3.7);
    EXPECT_EQ(data[48], 0.0);
    EXPECT_EQ(data[49], 3.5);
    EXPECT_EQ(data[50], 3.8);
    EXPECT_EQ(data[51], 0.0);
    EXPECT_EQ(data[52], 3.6);
    EXPECT_EQ(data[53], 3.9);

    // (2,1) starts at (2 + 1 * 4) * 9 = 54
    EXPECT_EQ(data[54], 2.1);
    EXPECT_EQ(data[55], 2.4);
    EXPECT_EQ(data[56], 2.7);
    EXPECT_EQ(data[57], 0.0);
    EXPECT_EQ(data[58], 2.5);
    EXPECT_EQ(data[59], 2.8);
    EXPECT_EQ(data[60], 0.0);
    EXPECT_EQ(data[61], 2.6);
    EXPECT_EQ(data[62], 2.9);

    // mutable access
    double *d = x->data_pointer();
    d[49]     = 4.1;

    // (1,1) starts at (1 + 1 * 4) * 9 = 45
    EXPECT_EQ(data[45], 3.1);
    EXPECT_EQ(data[46], 3.4);
    EXPECT_EQ(data[47], 3.7);
    EXPECT_EQ(data[48], 0.0);
    EXPECT_EQ(data[49], 4.1);
    EXPECT_EQ(data[50], 3.8);
    EXPECT_EQ(data[51], 0.0);
    EXPECT_EQ(data[52], 3.6);
    EXPECT_EQ(data[53], 3.9);
}

//---------------------------------------------------------------------------//

TEST_F(SDM_Face_Field_Test, Swap)
{
    if (nodes != 1) return;

    EXPECT_EQ(4, x->abscissa());
    EXPECT_EQ(2, x->ordinate());

    fill(x->begin_data(), x->end_data(), 1.1);

    SDM_Face_Field z(*mesh, X, 3);

    fill(z.begin_data(), z.end_data(), 2.1);

    x->swap(z);

    const double *zp = z.data_pointer();
    const double *xp = x->data_pointer();
    for (int i = 0; i < z.data_size(); ++i)
    {
        EXPECT_EQ(zp[i], 1.1);
        EXPECT_EQ(xp[i], 2.1);
    }
}

//---------------------------------------------------------------------------//

TEST_F(SDM_Face_Field_Test, Copy)
{
    if (nodes != 1) return;

    EXPECT_EQ(4, x->abscissa());
    EXPECT_EQ(2, x->ordinate());

    fill(x->begin_data(), x->end_data(), 1.1);

    SDM_Face_Field z(*mesh, X, 3);

    z.fast_copy(*x);

    fill(x->begin_data(), x->end_data(), 2.1);

    const double *zp = z.data_pointer();
    const double *xp = x->data_pointer();
    for (int i = 0; i < z.data_size(); ++i)
    {
        EXPECT_EQ(zp[i], 1.1);
        EXPECT_EQ(xp[i], 2.1);
    }

    vector<double> zd(72, 3.1);
    z.fast_copy(&zd[0], &zd[0] + zd.size());
    for (int i = 0; i < z.data_size(); ++i)
    {
        EXPECT_EQ(zp[i], 3.1);
        EXPECT_EQ(xp[i], 2.1);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstSDM_Face_Field.cc
//---------------------------------------------------------------------------//
