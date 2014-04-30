//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstFV_Gather.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 30 11:05:03 2012
 * \brief  FV_Gather test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "gtest/utils_gtest.hh"

#include "comm/SpinLock.hh"
#include "utils/Definitions.hh"
#include "mesh/Partitioner.hh"
#include "../Moment_Coefficients.hh"
#include "../Dimensions.hh"
#include "../SDM_Face_Field.hh"
#include "../FV_Gather.hh"
#include "Test_XS.hh"

using profugus::FV_Gather;
using namespace std;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class FV_Gather_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::Partitioner                         Partitioner;
    typedef Partitioner::RCP_ParameterList                RCP_ParameterList;
    typedef Partitioner::ParameterList                    ParameterList;
    typedef Partitioner::RCP_Mesh                         RCP_Mesh;
    typedef Partitioner::RCP_Indexer                      RCP_Indexer;
    typedef FV_Gather::RCP_Face_Field                     RCP_Face_Field;
    typedef profugus::Moment_Coefficients                 Moment_Coefficients_t;
    typedef FV_Gather::RCP_Moment_Coefficients            RCP_Moment_Coefficients;
    typedef profugus::Moment_Coefficients::Mat_DB_t       Mat_DB_t;
    typedef profugus::Moment_Coefficients::RCP_Mat_DB     RCP_Mat_DB;
    typedef profugus::Moment_Coefficients::RCP_Dimensions RCP_Dimensions;
    typedef FV_Gather::Face_Field_t::Serial_Matrix        Serial_Matrix;
    typedef Teuchos::RCP<Serial_Matrix>                   RCP_Serial_Matrix;

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
        db->set("num_cells_k", 3);

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
    }

    void problem(int Ng)
    {
        using def::I; using def::J; using def::K;

        // make 4 materials
        vector<double> f(5);
        vector<int>    ids(5);
        vector<int>    matids(mesh->num_cells(), 0);

        ids[0] = 100; f[0] = 0.8;
        ids[1] = 150; f[1] = 0.9;
        ids[2] = 200; f[2] = 1.1;
        ids[3] = 250; f[3] = 1.5;
        ids[4] = 1;   f[4] = 1.0;

        vector<int> gids(6*4*3, 0);

        for (int cell = 0; cell < 72; ++cell)
        {
            gids[cell] = 1;
        }

        for (int k = 0; k < 3; ++k)
        {
            for (int j = 0; j < 2; ++j)
            {
                int i      = 2;
                int cell   = i + j * 6 + k * 24;
                gids[cell] = 100;

                i          = 3;
                cell       = i + j * 6 + k * 24;
                gids[cell] = 150;

            }

            for (int j = 2; j < 4; ++j)
            {
                int i      = 2;
                int cell   = i + j * 6 + k * 24;
                gids[cell] = 200;

                i          = 3;
                cell       = i + j * 6 + k * 24;
                gids[cell] = 250;
            }
        }

        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_dim(I); ++i)
                {
                    int global    = indexer->l2g(i, j, k);
                    int local     = indexer->l2l(i, j, k);
                    matids[local] = gids[global];
                }
            }
        }

        if (Ng == 1)
            mat = one_grp::make_mat(3, ids, f, matids);
        if (Ng == 3)
            mat = three_grp::make_mat(3, ids, f, matids);

        EXPECT_FALSE(mat.is_null());

        RCP_ParameterList db = Teuchos::rcp(new ParameterList("MC"));
        dims                 = Teuchos::rcp(new profugus::Dimensions(7));
        mom_coeff            = Teuchos::rcp(new Moment_Coefficients_t(
                                                db, dims, mat));

        EXPECT_FALSE(dims.is_null());
        EXPECT_FALSE(mom_coeff.is_null());
    }

  protected:
    // >>> Data that get re-initialized between tests

    RCP_Mesh    mesh;
    RCP_Indexer indexer;
    RCP_Mat_DB  mat;

    RCP_Dimensions          dims;
    RCP_Moment_Coefficients mom_coeff;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(FV_Gather_Test, 1_PE_Test)
{
    if (nodes > 1) return;

    problem(3);

    FV_Gather g(mesh, mom_coeff, *indexer);
    g.gather(1);

    RCP_Face_Field lox = g.low_side_D(X);
    RCP_Face_Field loy = g.low_side_D(Y);
    RCP_Face_Field hix = g.high_side_D(X);
    RCP_Face_Field hiy = g.high_side_D(Y);

    EXPECT_TRUE(lox.is_null());
    EXPECT_TRUE(loy.is_null());
    EXPECT_TRUE(hix.is_null());
    EXPECT_TRUE(hiy.is_null());
}

//---------------------------------------------------------------------------//

TEST_F(FV_Gather_Test, 2_PE_1Grp_Test)
{
    if (nodes != 2) return;

    problem(1);

    FV_Gather g(mesh, mom_coeff, *indexer);
    g.gather(1);

    RCP_Face_Field lox = g.low_side_D(X);
    RCP_Face_Field loy = g.low_side_D(Y);
    RCP_Face_Field hix = g.high_side_D(X);
    RCP_Face_Field hiy = g.high_side_D(Y);

    if (node == 0)
    {
        EXPECT_TRUE(lox.is_null());
        EXPECT_TRUE(loy.is_null());
        EXPECT_FALSE(hix.is_null());
        EXPECT_TRUE(hiy.is_null());

        EXPECT_EQ(hix->abscissa(), 4);
        EXPECT_EQ(hix->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 7.0 / (0.9 * 0.97);
                Serial_Matrix m = hix->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }

            for (int a = 2; a < 4; ++a)
            {
                double ref      = 1.0 / 7.0 / (1.5 * 0.97);
                Serial_Matrix m = hix->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);

            }
        }
    }

    if (node == 1)
    {
        EXPECT_FALSE(lox.is_null());
        EXPECT_TRUE(loy.is_null());
        EXPECT_TRUE(hix.is_null());
        EXPECT_TRUE(hiy.is_null());

        EXPECT_EQ(lox->abscissa(), 4);
        EXPECT_EQ(lox->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 7.0 / (0.8 * 0.97);
                Serial_Matrix m = lox->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }

            for (int a = 2; a < 4; ++a)
            {
                double ref      = 1.0 / 7.0 / (1.1 * 0.97);
                Serial_Matrix m = lox->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);

            }
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(FV_Gather_Test, 4_PE_1Grp_Test)
{
    if (nodes != 4) return;

    problem(1);

    FV_Gather g(mesh, mom_coeff, *indexer);
    g.gather(2);

    RCP_Face_Field lox = g.low_side_D(X);
    RCP_Face_Field loy = g.low_side_D(Y);
    RCP_Face_Field hix = g.high_side_D(X);
    RCP_Face_Field hiy = g.high_side_D(Y);

    if (node == 0)
    {
        EXPECT_TRUE(lox.is_null());
        EXPECT_TRUE(loy.is_null());
        EXPECT_FALSE(hix.is_null());
        EXPECT_FALSE(hiy.is_null());

        EXPECT_EQ(hix->abscissa(), 2);
        EXPECT_EQ(hix->ordinate(), 3);

        EXPECT_EQ(hiy->abscissa(), 3);
        EXPECT_EQ(hiy->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 11.0 / (0.9);
                Serial_Matrix m = hix->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 11.0;
                Serial_Matrix m = hiy->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }

            {
                double ref      = 1.0 / 11.0 / (1.1);
                Serial_Matrix m = hiy->view(2, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }
    }

    if (node == 1)
    {
        EXPECT_FALSE(lox.is_null());
        EXPECT_TRUE(loy.is_null());
        EXPECT_TRUE(hix.is_null());
        EXPECT_FALSE(hiy.is_null());

        EXPECT_EQ(lox->abscissa(), 2);
        EXPECT_EQ(lox->ordinate(), 3);

        EXPECT_EQ(hiy->abscissa(), 3);
        EXPECT_EQ(hiy->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 11.0 / (0.8);
                Serial_Matrix m = lox->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }

        for (int o = 0; o < 3; ++o)
        {
            {
                double ref      = 1.0 / 11.0 / (1.5);
                Serial_Matrix m = hiy->view(0, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }

            for (int a = 1; a < 3; ++a)
            {
                double ref      = 1.0 / 11.0;
                Serial_Matrix m = hiy->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }
    }

    if (node == 2)
    {
        EXPECT_TRUE(lox.is_null());
        EXPECT_FALSE(loy.is_null());
        EXPECT_FALSE(hix.is_null());
        EXPECT_TRUE(hiy.is_null());

        EXPECT_EQ(hix->abscissa(), 2);
        EXPECT_EQ(hix->ordinate(), 3);

        EXPECT_EQ(loy->abscissa(), 3);
        EXPECT_EQ(loy->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 11.0 / (1.5);
                Serial_Matrix m = hix->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 11.0;
                Serial_Matrix m = loy->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }

            {
                double ref      = 1.0 / 11.0 / (0.8);
                Serial_Matrix m = loy->view(2, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }
    }

    if (node == 3)
    {
        EXPECT_FALSE(lox.is_null());
        EXPECT_FALSE(loy.is_null());
        EXPECT_TRUE(hix.is_null());
        EXPECT_TRUE(hiy.is_null());

        EXPECT_EQ(lox->abscissa(), 2);
        EXPECT_EQ(lox->ordinate(), 3);

        EXPECT_EQ(loy->abscissa(), 3);
        EXPECT_EQ(loy->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                double ref      = 1.0 / 11.0 / (1.1);
                Serial_Matrix m = lox->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }

        for (int o = 0; o < 3; ++o)
        {
            {
                double ref      = 1.0 / 11.0 / (0.9);
                Serial_Matrix m = loy->view(0, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }

            for (int a = 1; a < 3; ++a)
            {
                double ref      = 1.0 / 11.0;
                Serial_Matrix m = loy->view(a, o);

                EXPECT_NEAR(ref, m(0, 0), 1.0e-6);
            }
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(FV_Gather_Test, 2_PE_3Grp_Test)
{
    if (nodes != 2) return;

    problem(3);

    FV_Gather g(mesh, mom_coeff, *indexer);
    g.gather(0);

    RCP_Face_Field lox = g.low_side_D(X);
    RCP_Face_Field loy = g.low_side_D(Y);
    RCP_Face_Field hix = g.high_side_D(X);
    RCP_Face_Field hiy = g.high_side_D(Y);

    if (node == 0)
    {
        EXPECT_TRUE(lox.is_null());
        EXPECT_TRUE(loy.is_null());
        EXPECT_FALSE(hix.is_null());
        EXPECT_TRUE(hiy.is_null());

        EXPECT_EQ(hix->abscissa(), 4);
        EXPECT_EQ(hix->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                Serial_Matrix m = hix->view(a, o);
                EXPECT_NEAR(1.335948586285846, m(0, 0), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 1), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 2), 1.0e-6);
                EXPECT_NEAR(0.004161777972805, m(1, 0), 1.0e-6);
                EXPECT_NEAR(0.617328172558875, m(1, 1), 1.0e-6);
                EXPECT_NEAR(0.000032010454048, m(1, 2), 1.0e-6);
                EXPECT_NEAR(0.000052323586050, m(2, 0), 1.0e-6);
                EXPECT_NEAR(0.007761303935283, m(2, 1), 1.0e-6);
                EXPECT_NEAR(0.246217025507847, m(2, 2), 1.0e-6);
            }

            for (int a = 2; a < 4; ++a)
            {
                Serial_Matrix m = hix->view(a, o);
                EXPECT_NEAR(0.801569151771508, m(0, 0), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 1), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 2), 1.0e-6);
                EXPECT_NEAR(0.002497066783683, m(1, 0), 1.0e-6);
                EXPECT_NEAR(0.370396903535325, m(1, 1), 1.0e-6);
                EXPECT_NEAR(0.000019206272429, m(1, 2), 1.0e-6);
                EXPECT_NEAR(0.000031394151630, m(2, 0), 1.0e-6);
                EXPECT_NEAR(0.004656782361170, m(2, 1), 1.0e-6);
                EXPECT_NEAR(0.147730215304708, m(2, 2), 1.0e-6);
            }
        }
    }

    if (node == 1)
    {
        EXPECT_FALSE(lox.is_null());
        EXPECT_TRUE(loy.is_null());
        EXPECT_TRUE(hix.is_null());
        EXPECT_TRUE(hiy.is_null());

        EXPECT_EQ(lox->abscissa(), 4);
        EXPECT_EQ(lox->ordinate(), 3);

        for (int o = 0; o < 3; ++o)
        {
            for (int a = 0; a < 2; ++a)
            {
                Serial_Matrix m = lox->view(a, o);
                EXPECT_NEAR(1.502942159571577, m(0, 0), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 1), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 2), 1.0e-6);
                EXPECT_NEAR(0.004682000219405, m(1, 0), 1.0e-6);
                EXPECT_NEAR(0.694494194128734, m(1, 1), 1.0e-6);
                EXPECT_NEAR(0.000036011760804, m(1, 2), 1.0e-6);
                EXPECT_NEAR(0.000058864034306, m(2, 0), 1.0e-6);
                EXPECT_NEAR(0.008731466927193, m(2, 1), 1.0e-6);
                EXPECT_NEAR(0.276994153696328, m(2, 2), 1.0e-6);
            }

            for (int a = 2; a < 4; ++a)
            {
                Serial_Matrix m = lox->view(a, o);
                EXPECT_NEAR(1.093048843324783, m(0, 0), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 1), 1.0e-6);
                EXPECT_NEAR(0.0, m(0, 2), 1.0e-6);
                EXPECT_NEAR(0.003405091068658, m(1, 0), 1.0e-6);
                EXPECT_NEAR(0.505086686639079, m(1, 1), 1.0e-6);
                EXPECT_NEAR(0.000026190371494, m(1, 2), 1.0e-6);
                EXPECT_NEAR(0.000042810206768, m(2, 0), 1.0e-6);
                EXPECT_NEAR(0.006350157765231, m(2, 1), 1.0e-6);
                EXPECT_NEAR(0.201450293597330, m(2, 2), 1.0e-6);
            }
        }
    }
}

//---------------------------------------------------------------------------//
//                        end of tstFV_Gather.cc
//---------------------------------------------------------------------------//
