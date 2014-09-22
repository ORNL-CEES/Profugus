//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstLinear_System_FV.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 17 10:20:25 2014
 * \brief  Linear_System_FV unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <fstream>
#include <iomanip>

#include "gtest/utils_gtest.hh"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_RCP.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "utils/Definitions.hh"
#include "xs/Mat_DB.hh"
#include "mesh/Partitioner.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "../Isotropic_Source.hh"
#include "../Dimensions.hh"
#include "../Linear_System_FV.hh"
#include "../MatrixTraits.hh"
#include "../VectorTraits.hh"
#include "Test_XS.hh"

using namespace std;

int node, nodes;

typedef profugus::Isotropic_Source               External_Source;
typedef profugus::Mat_DB                         Mat_DB_t;
typedef Teuchos::RCP<Mat_DB_t>                   RCP_Mat_DB;
typedef Mat_DB_t::XS_t                           XS;
typedef Mat_DB_t::RCP_XS                         RCP_XS;
typedef profugus::Partitioner                    Partitioner;
typedef Partitioner::Array_Dbl                   Array_Dbl;
typedef Teuchos::RCP<profugus::Dimensions>       RCP_Dimensions;
typedef Teuchos::RCP<profugus::Mesh>             RCP_Mesh;
typedef Teuchos::RCP<profugus::Global_Mesh_Data> RCP_Global_Data;
typedef Teuchos::RCP<profugus::LG_Indexer>       RCP_Indexer;
typedef Teuchos::RCP<profugus::FV_Bnd_Indexer>   RCP_Bnd_Indexer;
typedef Teuchos::RCP<Teuchos::ParameterList>     RCP_ParameterList;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
template <class T>
class MatrixTest : public testing::Test
{
  protected:
    typedef profugus::Linear_System_FV<T>             Linear_System;
    typedef Teuchos::RCP<Linear_System>               RCP_Linear_System;
    typedef typename Linear_System::Matrix_t          Matrix_t;
    typedef typename Linear_System::Vector_t          Vector_t;

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
        db = Teuchos::rcp(new Partitioner::ParameterList("test"));

        if (cx.size() == 0)
        {
            db->set("delta_x", 1.0);
            db->set("delta_y", 1.0);
            db->set("delta_z", 1.0);

            db->set("num_cells_i", 4);
            db->set("num_cells_j", 4);
            db->set("num_cells_k", 4);
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

        Partitioner p(db);
        p.build();

        mesh    = p.get_mesh();
        indexer = p.get_indexer();
        data    = p.get_global_data();
    }

    void make_data()
    {
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
            dim    = Teuchos::rcp(new profugus::Dimensions(eqn_order));
            system = Teuchos::rcp(new Linear_System(
                                      db, dim, mat, mesh, indexer, data));
        }
        TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
    }

    void make_data(const vector<int>    &matids,
                   const vector<double> &f,
                   const vector<int>    &cell2mid)
    {
        if (num_groups == 1)
            mat = one_grp::make_mat(3, matids, f, cell2mid);
        else if (num_groups == 2)
            mat = two_grp::make_mat(3, matids, f, cell2mid);
        else
            mat = three_grp::make_mat(3, matids, f, cell2mid);

        EXPECT_FALSE(mat.is_null());
        EXPECT_FALSE(mesh.is_null());
        EXPECT_FALSE(indexer.is_null());
        EXPECT_FALSE(data.is_null());

        bool success = true, verbose = true;
        try
        {
            dim    = Teuchos::rcp(new profugus::Dimensions(eqn_order));
            system = Teuchos::rcp(new Linear_System(
                                      db, dim, mat, mesh, indexer, data));
        }
        TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
    }

  protected:

    RCP_Mesh        mesh;
    RCP_Indexer     indexer;
    RCP_Global_Data data;

    RCP_ParameterList db;
    RCP_Mat_DB        mat;
    RCP_Dimensions    dim;

    RCP_Linear_System system;

    int num_groups, eqn_order;

    Array_Dbl cx, cy, cz;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(MatrixTest, MyTypes);

TYPED_TEST(MatrixTest, SP7_3Grp_Refl_Graph)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;

    this->build(7, 3);
    this->make_data();
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Dimensions dim = this->dim;
    RCP_Global_Data data = this->data;
    RCP_Linear_System system = this->system;
    EXPECT_EQ(4, dim->num_equations());

    // make the matrix
    system->build_Matrix();

    // get the graph
    Teuchos::RCP<const Matrix_t> A = Teuchos::rcp_dynamic_cast<const Matrix_t>(
        system->get_Operator());

    EXPECT_EQ(4 * 3 * mesh->num_cells(), MatrixTraits::local_rows(A));
    EXPECT_EQ(4 * 3 * data->num_cells(), MatrixTraits::global_rows(A));
    EXPECT_EQ(4 * 3 * data->num_cells(), MatrixTraits::global_columns(A));
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP5_3Grp_Refl_Matrix)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;

    this->build(5, 3);
    this->make_data();
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Dimensions dim = this->dim;
    RCP_Global_Data data = this->data;
    RCP_Linear_System system = this->system;
    EXPECT_EQ(3, dim->num_equations());

    // make the matrix
    system->build_Matrix();
    Teuchos::RCP<const Matrix_t> A = Teuchos::rcp_dynamic_cast<const Matrix_t>(
        system->get_Operator());

    EXPECT_EQ(3 * 3 * mesh->num_cells(), MatrixTraits::local_rows(A));
    EXPECT_EQ(3 * 3 * data->num_cells(), MatrixTraits::global_rows(A));
    EXPECT_EQ(3 * 3 * data->num_cells(), MatrixTraits::global_columns(A));
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP1_2Grp_Refl_Matrix)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;
    typedef typename TypeParam::MV                  MV;
    typedef typename TypeParam::OP                  OP;
    typedef Anasazi::OperatorTraits<double,MV,OP>   OPT;

    using def::I; using def::J; using def::K;
    this->build(1, 2);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    // 2 materials
    vector<int>    ids(2, 0);
    vector<double> f(2, 0.0);
    vector<int>    matids(mesh->num_cells(), 0);
    ids[0] = 9;  f[0] = 0.9;
    ids[1] = 11; f[1] = 1.1;

    vector<int> gids(4*4*4, 0);
    for (int k = 0; k < 4; ++k)
    {
        for (int j = 0; j < 4; ++j)
        {
            gids[indexer->g2g(0, j, k)] = 9;
            gids[indexer->g2g(1, j, k)] = 11;
            gids[indexer->g2g(2, j, k)] = 9;
            gids[indexer->g2g(3, j, k)] = 11;
        }
    }

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                matids[indexer->l2l(i, j, k)] = gids[indexer->l2g(i, j, k)];
            }
        }
    }

    this->make_data(ids, f, matids);
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    int num_groups = this->num_groups;

    EXPECT_EQ(1, dim->num_equations());

    // make the matrix
    system->build_Matrix();
    Teuchos::RCP<const Matrix_t> A = Teuchos::rcp_dynamic_cast<const Matrix_t>(
        system->get_Operator());

    EXPECT_EQ(2 * mesh->num_cells(), MatrixTraits::local_rows(A));
    EXPECT_EQ(2 * data->num_cells(), MatrixTraits::global_rows(A));
    EXPECT_EQ(2 * data->num_cells(), MatrixTraits::global_columns(A));

    // make a vector
    Teuchos::RCP<Vector_t> x = VectorTraits::build_vector(system->get_Map());
    Teuchos::RCP<Vector_t> y = VectorTraits::build_vector(system->get_Map());
    EXPECT_EQ(mesh->num_cells() * 2, VectorTraits::local_length(x));
    EXPECT_EQ(mesh->num_cells() * 2, VectorTraits::local_length(y));
    Teuchos::ArrayView<double> x_data = VectorTraits::get_data_nonconst(x,0);
    Teuchos::ArrayView<double> y_data = VectorTraits::get_data_nonconst(y,0);

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int global = indexer->l2g(i, j, k);
                int local  = indexer->l2l(i, j, k);

                double gv = static_cast<double>(global);

                x_data[0 + local * num_groups] = gv + 0.1;
                x_data[1 + local * num_groups] = gv + 0.2;
            }
        }
    }

    // extract some blocks
    if (node == 0)
    {
        Teuchos::ArrayView<const int>    indices;
        Teuchos::ArrayView<const double> values;
        MatrixTraits::get_local_row_view(A,2,indices,values);
        EXPECT_EQ(5, indices.size());
        EXPECT_EQ(5, values.size());

        if (nodes == 1)
        {
            EXPECT_EQ(0,  MatrixTraits::global_col_id(A,indices[0]));
            EXPECT_EQ(2,  MatrixTraits::global_col_id(A,indices[1]));
            EXPECT_EQ(4,  MatrixTraits::global_col_id(A,indices[2]));
            EXPECT_EQ(10, MatrixTraits::global_col_id(A,indices[3]));
            EXPECT_EQ(34, MatrixTraits::global_col_id(A,indices[4]));

            EXPECT_SOFTEQ(-0.24691358, values[0], 1.0e-6);
            EXPECT_SOFTEQ( 1.49276094, values[1], 1.0e-6);
            EXPECT_SOFTEQ(-0.24691358, values[2], 1.0e-6);
            EXPECT_SOFTEQ(-0.22446689, values[3], 1.0e-6);
            EXPECT_SOFTEQ(-0.22446689, values[4], 1.0e-6);
        }
    }

    OPT::Apply(*A,*x,*y);

    double v[128] = {
        -5.688882,   -5.007769,   -3.884338,   -4.254782,   -4.541968,
        -5.332511,   -2.537424,   -4.699524,   -2.791488,   -5.131067,
        -0.786470,   -4.792025,   -1.644575,   -5.455809,    0.560443,
        -5.236767,   -0.991488,   -6.211067,    1.413530,   -6.112025,
        0.155425,   -6.535809,    2.760443,   -6.556767,    1.905905,
        -6.334365,    4.511397,   -6.649269,    3.052819,   -6.659107,
        5.858311,   -7.094011,    5.900693,   -5.500960,    8.507132,
        -6.403756,    7.047606,   -5.825702,    9.854046,   -6.848498,
        8.798086,   -5.624258,   11.605000,   -6.941000,    9.945000,
        -5.949000,   12.951914,   -7.385742,   10.598086,   -6.704258,
        13.805000,   -8.261000,   11.745000,   -7.029000,   15.151914,
        -8.705742,   13.495480,   -6.827556,   16.902868,   -8.798244,
        14.642394,   -7.152298,   18.249781,   -9.242986,   13.100693,
        -9.820960,   17.307132,  -11.683756,   14.247606,  -10.145702,
        18.654046,  -12.128498,   15.998086,   -9.944258,   20.405000,
        -12.221000,   17.145000,  -10.269000,   21.751914,  -12.665742,
        17.798086,  -11.024258,   22.605000,  -13.541000,   18.945000,
        -11.349000,   23.951914,  -13.985742,   20.695480,  -11.147556,
        25.702868,  -14.078244,   21.842394,  -11.472298,   27.049781,
        -14.522986,   24.690267,  -10.314151,   29.698603,  -13.832731,
        25.837181,  -10.638893,   31.045516,  -14.277473,   27.587661,
        -10.437449,   32.796470,  -14.369975,   28.734575,  -10.762191,
        34.143384,  -14.814717,   29.387661,  -11.517449,   34.996470,
        -15.689975,   30.534575,  -11.842191,   36.343384,  -16.134717,
        32.285055,  -11.640747,   38.094338,  -16.227218,   33.431968,
        -11.965489,   39.441251,  -16.671960};

    double eps = 1.0e-5;
    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int g = 0; g < num_groups; ++g)
                {
                    int local  = g + indexer->l2l(i, j, k) * 2;
                    int global = g + indexer->l2g(i, j, k) * 2;
                    EXPECT_SOFTEQ(v[global], y_data[local], eps);
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP3_2Grp_Refl_Matrix)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;
    typedef typename TypeParam::MV                  MV;
    typedef typename TypeParam::OP                  OP;
    typedef Anasazi::OperatorTraits<double,MV,OP>   OPT;

    using def::I; using def::J; using def::K;

    Array_Dbl &cx = this->cx;
    Array_Dbl &cy = this->cy;
    Array_Dbl &cz = this->cz;

    // make non-uniform mesh
    cx.resize(5);
    cy.resize(5);
    cz.resize(5);

    cx[0] = 0.0;
    cx[1] = 0.8;
    cx[2] = 1.7;
    cx[3] = 2.7;
    cx[4] = 3.8;

    cy[0] = 0.0;
    cy[1] = 0.7;
    cy[2] = 1.5;
    cy[3] = 2.4;
    cy[4] = 3.4;

    cz[0] = 0.0;
    cz[1] = 0.6;
    cz[2] = 1.3;
    cz[3] = 2.1;
    cz[4] = 3.0;

    // build the mesh and data
    this->build(3, 2);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    // 2 materials
    vector<int>    ids(2, 0);
    vector<double> f(2, 0.0);
    vector<int>    matids(mesh->num_cells(), 0);
    ids[0] = 9;  f[0] = 0.9;
    ids[1] = 11; f[1] = 1.1;

    vector<int> gids(4*4*4, 0);
    for (int k = 0; k < 4; ++k)
    {
        for (int j = 0; j < 4; ++j)
        {
            gids[indexer->g2g(0, j, k)] = 9;
            gids[indexer->g2g(1, j, k)] = 11;
            gids[indexer->g2g(2, j, k)] = 9;
            gids[indexer->g2g(3, j, k)] = 11;
        }
    }

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                matids[indexer->l2l(i, j, k)] = gids[indexer->l2g(i, j, k)];
            }
        }
    }

    this->make_data(ids, f, matids);
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    int num_groups = this->num_groups;

    EXPECT_EQ(2, dim->num_equations());
    EXPECT_EQ("reflect", db->get<string>("boundary"));

    // make the matrix
    this->system->build_Matrix();
    Teuchos::RCP<const Matrix_t> A = Teuchos::rcp_dynamic_cast<const Matrix_t>(
        system->get_Operator());

    EXPECT_EQ(2 * 2 * mesh->num_cells(), MatrixTraits::local_rows(A));
    EXPECT_EQ(2 * 2 * data->num_cells(), MatrixTraits::global_rows(A));
    EXPECT_EQ(2 * 2 * data->num_cells(), MatrixTraits::global_columns(A));

    // make a vector
    Teuchos::RCP<Vector_t> x = VectorTraits::build_vector(system->get_Map());
    Teuchos::RCP<Vector_t> y = VectorTraits::build_vector(system->get_Map());
    EXPECT_EQ(mesh->num_cells() * 2 * 2, VectorTraits::local_length(x));
    EXPECT_EQ(mesh->num_cells() * 2 * 2, VectorTraits::local_length(y));
    Teuchos::ArrayView<double> x_data = VectorTraits::get_data_nonconst(x,0);
    Teuchos::ArrayView<double> y_data = VectorTraits::get_data_nonconst(y,0);

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int global = indexer->l2g(i, j, k);
                int local  = indexer->l2l(i, j, k);

                for (int n = 0; n < 2; ++n)
                {
                    double gv = static_cast<double>(global + n);
                    x_data[0 + n * num_groups + local * num_groups * 2] = gv + 0.1;
                    x_data[1 + n * num_groups + local * num_groups * 2] = gv + 0.2;
                }
            }
        }
    }

    // extract some blocks
    if (node == 0)
    {
        Teuchos::ArrayView<const int>    indices;
        Teuchos::ArrayView<const double> values;
        MatrixTraits::get_local_row_view(A,2,indices,values);
        EXPECT_EQ(5, indices.size());
        EXPECT_EQ(5, values.size());

        if (nodes == 1)
        {
            EXPECT_EQ(0,  MatrixTraits::global_col_id(A,indices[0]));
            EXPECT_EQ(2,  MatrixTraits::global_col_id(A,indices[1]));
            EXPECT_EQ(6,  MatrixTraits::global_col_id(A,indices[2]));
            EXPECT_EQ(18, MatrixTraits::global_col_id(A,indices[3]));
            EXPECT_EQ(66, MatrixTraits::global_col_id(A,indices[4]));

            EXPECT_SOFTEQ(-0.30000000, values[0], 1.0e-6);
            EXPECT_SOFTEQ( 1.52962520, values[1], 1.0e-6);
            EXPECT_SOFTEQ(-0.14207855, values[2], 1.0e-6);
            EXPECT_SOFTEQ(-0.20567562, values[3], 1.0e-6);
            EXPECT_SOFTEQ(-0.27687103, values[4], 1.0e-6);
        }
    }

    OPT::Apply(*A,*x,*y);

    double v[256] = {
          -13.9915783590,   -11.7723247232,    -4.4292174299,    -3.8797514745,
          -11.0535624733,    -9.3862595921,    -2.3660858131,    -1.6364723952,
          -13.3033506889,   -11.6138698313,    -3.0664164519,    -1.9853602529,
          -10.5046905378,    -9.4474139732,    -0.8154834802,     0.5777119859,
          -11.0861294304,   -10.1224461701,    -1.1018249945,     0.4383165711,
           -8.4339527439,    -8.1818135032,     1.3341039573,     3.3187145916,
          -10.3979017604,    -9.9639912783,     0.2609759835,     2.3327077927,
           -7.8850808084,    -8.2429678843,     2.8847062902,     5.5328989727,
          -10.5503045000,   -10.5383936667,     1.2929165962,     3.9361013437,
           -7.7531262857,    -8.6675887277,     4.2712157638,     7.6027607187,
           -9.8620768299,   -10.3799387748,     2.6557175742,     5.8304925653,
           -7.2042543502,    -8.7287431087,     5.8218180967,     9.8169450998,
           -8.9461534111,   -10.0229798969,     4.1081364125,     7.8037043146,
           -6.1982147887,    -8.3913410978,     7.5523552095,    12.1893853717,
           -8.2579257411,    -9.8645250050,     5.4709373905,     9.6980955363,
           -5.6493428533,    -8.4524954789,     9.1029575424,    14.4035697528,
            0.9500638123,    -2.2785853944,    10.1869974188,    14.5417386937,
            2.1411144547,    -2.2004728685,    13.5037465379,    19.1245448131,
            1.6382914824,    -2.1201305025,    11.5497983969,    16.4361299153,
            2.6899863902,    -2.2616272496,    15.0543488708,    21.3387291942,
            3.8555127409,    -0.6287068413,    13.5143898543,    18.8598067393,
            4.7607241841,    -0.9960267796,    17.2039363083,    24.0797317999,
            4.5437404109,    -0.4702519495,    14.8771908323,    20.7541979609,
            5.3095961196,    -1.0571811607,    18.7545386412,    26.2939161811,
            4.3913376713,    -1.0446543379,    15.9091314450,    22.3575915118,
            5.4415506424,    -1.4818020041,    20.1410481148,    28.3637779270,
            5.0795653413,    -0.8861994460,    17.2719324230,    24.2519827335,
            5.9904225778,    -1.5429563852,    21.6916504477,    30.5779623081,
            5.9954887602,    -0.5292405681,    18.7243512613,    26.2251944828,
            6.9964621393,    -1.2055543743,    23.4221875605,    32.9504025800,
            6.6837164302,    -0.3707856762,    20.0871522393,    28.1195857044,
            7.5453340748,    -1.2667087553,    24.9727898934,    35.1645869611,
            2.9244424717,    -4.0896398965,    19.6994788154,    28.4744030545,
            4.7262121457,    -4.2640629157,    25.1977969734,    36.2128863608,
            3.6126701417,    -3.9311850047,    21.0622797934,    30.3687942761,
            5.2750840811,    -4.3252172968,    26.7483993063,    38.4270707420,
            5.8298914002,    -2.4397613434,    23.0268712508,    32.7924711001,
            7.3458218751,    -3.0596168268,    28.8979867439,    41.1680733477,
            6.5181190703,    -2.2813064516,    24.3896722288,    34.6868623217,
            7.8946938106,    -3.1207712079,    30.4485890767,    43.3822577288,
            6.3657163306,    -2.8557088400,    25.4216128415,    36.2902558726,
            8.0266483333,    -3.5453920513,    31.8350985504,    45.4521194747,
            7.0539440007,    -2.6972539481,    26.7844138195,    38.1846470943,
            8.5755202688,    -3.6065464324,    33.3857008832,    47.6663038558,
            7.9698674195,    -2.3402950702,    28.2368326578,    40.1578588436,
            9.5815598303,    -3.2691444215,    35.1162379960,    50.0387441277,
            8.6580950896,    -2.1818401783,    29.5996336358,    42.0522500652,
           10.1304317657,    -3.3302988025,    36.6668403289,    52.2529285088,
           10.2017477605,    -1.2776301575,    31.2991179232,    44.2427603421,
           11.6500679881,    -2.5451458565,    38.5995219001,    54.8031584851,
           10.8899754305,    -1.1191752657,    32.6619189012,    46.1371515637,
           12.1989399235,    -2.6063002376,    40.1501242329,    57.0173428662,
           13.1071966890,     0.3722483956,    34.6265103587,    48.5608283877,
           14.2696777175,    -1.3406997677,    42.2997116705,    59.7583454719,
           13.7954243591,     0.5307032874,    35.9893113367,    50.4552196094,
           14.8185496529,    -1.4018541487,    43.8503140033,    61.9725298530,
           13.6430216194,    -0.0436991010,    37.0212519493,    52.0586131603,
           14.9505041757,    -1.8264749921,    45.2368234770,    64.0423915989,
           14.3312492895,     0.1147557908,    38.3840529274,    53.9530043819,
           15.4993761111,    -1.8876293732,    46.7874258099,    66.2565759800,
           15.2471727083,     0.4717146688,    39.8364717657,    55.9262161312,
           16.5054156727,    -1.5502273623,    48.5179629227,    68.6290162519,
           15.9354003784,     0.6301695606,    41.1992727437,    57.8206073529,
           17.0542876081,    -1.6113817434,    50.0685652555,    70.8432006331};

    double eps = 1.0e-6;
    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 2; ++n)
                {
                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + indexer->l2l(i, j, k) * 4;
                        int global = g + n * 2 + indexer->l2g(i, j, k) * 4;
                        EXPECT_EQ(local, system->index(
                                      g, n, indexer->l2l(i, j, k)));
                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP3_2Grp_Vac_Graph)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;

    using def::I; using def::J; using def::K;

    // build the mesh and data
    this->build(3, 2);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    // make vacuum boundary conditions
    db->set("boundary", string("vacuum"));

    // make the graph and data
    this->make_data();
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    EXPECT_EQ(2, dim->num_equations());

    // make the matrix
    system->build_Matrix();

    // get the graph
    Teuchos::RCP<const Matrix_t> A = Teuchos::rcp_dynamic_cast<const Matrix_t>(
        system->get_Operator());

    // number of faces
    int nf = 0;
    if (nodes == 1)
    {
        nf = 2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(J)) +
             2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(K)) +
             2 * (mesh->num_cells_dim(J) * mesh->num_cells_dim(K));
    }
    if (nodes == 2)
    {
        nf = 2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(J)) +
             2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(K)) +
             1 * (mesh->num_cells_dim(J) * mesh->num_cells_dim(K));
    }
    if (nodes == 4)
    {
        nf = 2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(J)) +
             1 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(K)) +
             1 * (mesh->num_cells_dim(J) * mesh->num_cells_dim(K));
    }

    EXPECT_EQ(2 * 2 * (mesh->num_cells() + nf), MatrixTraits::local_rows(A));
    EXPECT_EQ(2 * 2 * (data->num_cells() + 96), MatrixTraits::global_rows(A));
    EXPECT_EQ(2 * 2 * (data->num_cells() + 96), MatrixTraits::global_columns(A));
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP3_2Grp_Vac_Matrix)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;
    typedef typename TypeParam::MV                  MV;
    typedef typename TypeParam::OP                  OP;
    typedef Anasazi::OperatorTraits<double,MV,OP>   OPT;

    Array_Dbl &cx = this->cx;
    Array_Dbl &cy = this->cy;
    Array_Dbl &cz = this->cz;

    using def::I; using def::J; using def::K;

    // make non-uniform mesh
    cx.resize(5);
    cy.resize(5);
    cz.resize(5);

    cx[0] = 0.0;
    cx[1] = 0.8;
    cx[2] = 1.7;
    cx[3] = 2.7;
    cx[4] = 3.8;

    cy[0] = 0.0;
    cy[1] = 0.7;
    cy[2] = 1.5;
    cy[3] = 2.4;
    cy[4] = 3.4;

    cz[0] = 0.0;
    cz[1] = 0.6;
    cz[2] = 1.3;
    cz[3] = 2.1;
    cz[4] = 3.0;

    // build the mesh and data
    this->build(3, 2);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    // 2 materials
    vector<int>    ids(2, 0);
    vector<double> f(2, 0.0);
    vector<int>    matids(mesh->num_cells(), 0);
    ids[0] = 9;  f[0] = 0.9;
    ids[1] = 11; f[1] = 1.1;

    vector<int> gids(4*4*4, 0);
    for (int k = 0; k < 4; ++k)
    {
        for (int j = 0; j < 4; ++j)
        {
            gids[indexer->g2g(0, j, k)] = 9;
            gids[indexer->g2g(1, j, k)] = 11;
            gids[indexer->g2g(2, j, k)] = 9;
            gids[indexer->g2g(3, j, k)] = 11;
        }
    }

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                matids[indexer->l2l(i, j, k)] = gids[indexer->l2g(i, j, k)];
            }
        }
    }

    // make vacuum boundary conditions
    db->set("boundary", string("vacuum"));

    // make the graph and data
    this->make_data(ids, f, matids);
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    int num_groups = this->num_groups;
    EXPECT_EQ(2, dim->num_equations());

    // make the matrix
    system->build_Matrix();
    Teuchos::RCP<const Matrix_t> A = Teuchos::rcp_dynamic_cast<const Matrix_t>(
        system->get_Operator());

    // number of faces
    int nf = 0;
    if (nodes == 1)
    {
        nf = 2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(J)) +
             2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(K)) +
             2 * (mesh->num_cells_dim(J) * mesh->num_cells_dim(K));
    }
    if (nodes == 2)
    {
        nf = 2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(J)) +
             2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(K)) +
             1 * (mesh->num_cells_dim(J) * mesh->num_cells_dim(K));
    }
    if (nodes == 4)
    {
        nf = 2 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(J)) +
             1 * (mesh->num_cells_dim(I) * mesh->num_cells_dim(K)) +
             1 * (mesh->num_cells_dim(J) * mesh->num_cells_dim(K));
    }

    EXPECT_EQ(nf*2*2, system->bnd_unknowns());

    EXPECT_EQ(2 * 2 * (mesh->num_cells() + nf), MatrixTraits::local_rows(A));
    EXPECT_EQ(2 * 2 * (data->num_cells() + 96), MatrixTraits::global_rows(A));
    EXPECT_EQ(2 * 2 * (data->num_cells() + 96), MatrixTraits::global_columns(A));

    // make a vector
    Teuchos::RCP<Vector_t> x = VectorTraits::build_vector(system->get_Map());
    Teuchos::RCP<Vector_t> y = VectorTraits::build_vector(system->get_Map());
    EXPECT_EQ((mesh->num_cells() + nf) * 2 * 2, VectorTraits::local_length(x));
    EXPECT_EQ((mesh->num_cells() + nf) * 2 * 2, VectorTraits::local_length(y));
    Teuchos::ArrayView<double> x_data = VectorTraits::get_data_nonconst(x,0);
    Teuchos::ArrayView<double> y_data = VectorTraits::get_data_nonconst(y,0);

    VectorTraits::put_scalar(x,0.4);

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int global = indexer->l2g(i, j, k);
                int local  = indexer->l2l(i, j, k);

                for (int n = 0; n < 2; ++n)
                {
                    double gv = static_cast<double>(global + n);
                    x_data[0 + n * num_groups + local * num_groups * 2] = gv + 0.1;
                    x_data[1 + n * num_groups + local * num_groups * 2] = gv + 0.2;
                }
            }
        }
    }

    // extract some blocks
    if (node == 0)
    {
        Teuchos::ArrayView<const int>    indices;
        Teuchos::ArrayView<const double> values;
        MatrixTraits::get_local_row_view(A,2,indices,values);
        EXPECT_EQ(8, indices.size());
        EXPECT_EQ(8, values.size());

        if (nodes == 1)
        {
            EXPECT_EQ(0,   MatrixTraits::global_col_id(A,indices[0]));
            EXPECT_EQ(2,   MatrixTraits::global_col_id(A,indices[1]));
            EXPECT_EQ(6,   MatrixTraits::global_col_id(A,indices[2]));
            EXPECT_EQ(18,  MatrixTraits::global_col_id(A,indices[3]));
            EXPECT_EQ(66,  MatrixTraits::global_col_id(A,indices[4]));
            EXPECT_EQ(258, MatrixTraits::global_col_id(A,indices[5]));
            EXPECT_EQ(386, MatrixTraits::global_col_id(A,indices[6]));
            EXPECT_EQ(514, MatrixTraits::global_col_id(A,indices[7]));

            EXPECT_SOFTEQ(-0.30000000, values[0], 1.0e-6);
            EXPECT_SOFTEQ( 2.90768245, values[1], 1.0e-6);
            EXPECT_SOFTEQ(-0.14207855, values[2], 1.0e-6);
            EXPECT_SOFTEQ(-0.20567562, values[3], 1.0e-6);
            EXPECT_SOFTEQ(-0.27687103, values[4], 1.0e-6);
            EXPECT_SOFTEQ(-0.33743656, values[5], 1.0e-6);
            EXPECT_SOFTEQ(-0.44073347, values[6], 1.0e-6);
            EXPECT_SOFTEQ(-0.59988722, values[7], 1.0e-6);
        }
    }

    OPT::Apply(*A,*x,*y);

    double v[640] = {
          -15.0419642197,   -12.3850498086,    -3.4645773537,    -2.9093003135,
           -9.5393016998,    -7.8789322470,    -0.9186770354,    -0.2880536579,
           -8.8086401389,    -7.4666010848,    -0.2567405893,     0.5779517958,
           -3.6622173435,    -3.2628708937,     2.8750632864,     3.9119425310,
           -2.2745915551,    -2.2345015752,     3.3035967862,     4.3959686169,
           -2.5728728087,    -2.9642372980,     4.1317598162,     5.8227673784,
           -1.7102017146,    -2.2582006507,     4.2802203656,     5.9208284068,
            2.9559073818,     1.3481742699,     7.7884124926,     9.9021925058,
            7.7872202676,     5.6542574497,     9.4476335094,    11.1913262168,
            3.0961067858,     0.8986338015,     9.0321388921,    11.8335425496,
            4.9222548268,     2.6408640265,     9.0745108411,    11.5290597889,
           10.1089655954,     6.5048574988,    13.2729041445,    16.4267052438,
           25.3371113351,    20.1185689407,    18.7548328336,    20.7878666054,
           15.3406304539,    10.5328805424,    16.6972541503,    20.2916001619,
           20.1401843031,    15.0721833610,    17.4638827168,    20.3183543493,
           24.7354354470,    18.2156727075,    21.8755200872,    25.7093116509,
           31.9909781204,    24.9538221566,    23.1824369293,    26.0403728381,
           17.4415107075,    11.2176199146,    19.8863684974,    24.7701088355,
           21.4585652709,    15.2559874027,    19.7915142604,    23.7238967711,
           26.7608329718,    18.8346066238,    25.0349254885,    30.1616907698,
           20.7450875008,    14.1697306037,    20.4993266862,    25.0330663472,
            4.7607241841,    -0.9960267796,    17.2039363083,    24.0797317999,
            4.5437404109,    -0.4702519495,    14.8771908323,    20.7541979609,
           13.7317422828,     6.3173101440,    22.2154056541,    29.3507426726,
           24.7102677124,    16.7434774549,    24.2438145245,    29.7179773376,
            5.4415506424,    -1.4818020041,    20.1410481148,    28.3637779270,
            5.0795653413,    -0.8861994460,    17.2719324230,    24.2519827335,
           15.8966473600,     7.1253470490,    25.7366300366,    34.1485263665,
           44.9426766888,    33.5463943023,    34.6068153580,    40.2431142344,
           19.8808616904,    10.0658491194,    28.6700009716,    37.5815963186,
           22.9800127265,    13.8837252449,    26.7171058080,    33.9698738021,
           32.7179045917,    20.7495666402,    35.2030835779,    44.1908935530,
           65.5994095779,    50.7211726579,    45.1456388149,    50.9236708820,
           34.6856706766,    21.9337251868,    37.3500207042,    46.9329194109,
           41.3495756035,    29.0645605131,    36.3557311551,    43.8586899436,
           49.9413074172,    34.7259604137,    44.8350279994,    54.3794516127,
           36.4368872850,    24.3174534925,    35.4107930735,    43.7142355795,
            7.3458218751,    -3.0596168268,    28.8979867439,    41.1680733477,
            6.5181190703,    -2.2813064516,    24.3896722288,    34.6868623217,
           22.2531544500,     9.4289686146,    36.2459063937,    48.4940344878,
           40.4020674966,    26.8912003437,    39.1552809119,    48.3991465699,
            8.0266483333,    -3.5453920513,    31.8350985504,    45.4521194747,
            7.0539440007,    -2.6972539481,    26.7844138195,    38.1846470943,
           24.4180595273,    10.2370055196,    39.7671307762,    53.2918181817,
           69.4136259929,    51.3477347212,    52.9736321395,    61.9633265844,
           29.6488998976,    14.2643097787,    43.1911562750,    57.1564276900,
           33.7335409057,    19.7262882730,    39.6849375987,    50.9415812806,
           48.4222572752,    30.1232758173,    52.0606891854,    65.8206751917,
          136.8229152032,   109.3397733789,    82.1802240440,    89.0861329702,
           83.2600396568,    60.0113090182,    67.3629455453,    80.1523593507,
          100.2104671287,    76.9054976596,    68.5245644790,    77.7416352958,
          105.5604693484,    78.9451586126,    77.6206847600,    90.0375593852,
           92.4533752865,    69.6786542337,    66.4601089755,    76.6125451268,
           43.4780855267,    24.1710033847,    54.0138730445,    70.0804793904,
           50.1719929716,    32.3022307151,    50.5732362255,    63.3056551331,
           65.4302110526,    42.8012031936,    64.1345246138,    79.8451083634,
           99.1281695474,    74.6146287176,    71.2710629849,    82.2354323881,
           46.3758689344,    25.6179598597,    57.8235480819,    75.1319606482,
           53.4174319514,    34.2485108513,    54.0344439872,    67.7414161765,
           69.8120730794,    45.5419717981,    68.5283122272,    85.4103271880,
          139.6284916130,   109.0870082581,    89.6112307776,    99.7766317913,
           77.3980179644,    51.6224440957,    72.8792739053,    90.0901938178,
           91.5857924257,    66.6878982353,    71.4567843315,    84.4753697515,
          103.2161682930,    73.6230245018,    84.5215387351,   101.1931091524,
            0.3557613169,     0.2700274348,    -0.1222978080,    -0.1234361482,
           -2.3877229081,    -2.1217280433,    -1.2020948062,    -1.0731371225,
           -5.1312071331,    -4.5134835215,    -2.2818918043,    -2.0228380968,
           -7.8746913580,    -6.9052389997,    -3.3616888025,    -2.9725390711,
          -10.6181755830,    -9.2969944779,    -4.4414858007,    -3.9222400454,
          -13.3616598080,   -11.6887499560,    -5.5212827988,    -4.8719410197,
          -16.1051440329,   -14.0805054342,    -6.6010797970,    -5.8216419940,
          -18.8486282579,   -16.4722609124,    -7.6808767952,    -6.7713429682,
          -21.5921124829,   -18.8640163906,    -8.7606737933,    -7.7210439425,
          -24.3355967078,   -21.2557718687,    -9.8404707915,    -8.6707449168,
          -27.0790809328,   -23.6475273469,   -10.9202677897,    -9.6204458911,
          -29.8225651578,   -26.0392828251,   -12.0000647878,   -10.5701468654,
          -32.5660493827,   -28.4310383033,   -13.0798617860,   -11.5198478397,
          -35.3095336077,   -30.8227937814,   -14.1596587841,   -12.4695488140,
          -38.0530178326,   -33.2145492596,   -15.2394557823,   -13.4192497883,
          -40.7965020576,   -35.6063047378,   -16.3192527805,   -14.3689507626,
           -0.9519283747,    -0.8459737232,    -0.5276678794,    -0.4702858563,
           -2.5844148556,    -2.2691670656,    -1.1701917130,    -1.0353971799,
           -4.2169013366,    -3.6923604080,    -1.8127155467,    -1.6005085034,
           -5.8493878176,    -5.1155537504,    -2.4552393803,    -2.1656198270,
           -7.4818742985,    -6.5387470928,    -3.0977632139,    -2.7307311505,
           -9.1143607795,    -7.9619404352,    -3.7402870475,    -3.2958424741,
          -10.7468472605,    -9.3851337776,    -4.3828108811,    -3.8609537976,
          -12.3793337415,   -10.8083271199,    -5.0253347148,    -4.4260651211,
          -14.0118202224,   -12.2315204623,    -5.6678585484,    -4.9911764447,
          -15.6443067034,   -13.6547138047,    -6.3103823820,    -5.5562877682,
          -17.2767931844,   -15.0779071471,    -6.9529062156,    -6.1213990918,
          -18.9092796653,   -16.5011004895,    -7.5954300492,    -6.6865104153,
          -20.5417661463,   -17.9242938319,    -8.2379538829,    -7.2516217389,
          -22.1742526273,   -19.3474871743,    -8.8804777165,    -7.8167330624,
          -23.8067391083,   -20.7706805166,    -9.5230015501,    -8.3818443860,
          -25.4392255892,   -22.1938738590,   -10.1655253837,    -8.9469557095,
            0.3851557907,     0.2871742112,    -0.1492927330,    -0.1505936932,
           -0.2989337823,    -0.2968782247,    -0.3624474651,    -0.3330998296,
           -1.1825494807,    -1.0795432049,    -0.7663195891,    -0.6932799642,
           -1.5816017316,    -1.4151015651,    -0.8672876201,    -0.7771158695,
          -12.1564863806,   -10.6465651176,    -5.0855075817,    -4.4920838614,
          -10.5602773769,    -9.2426649482,    -4.4011687050,    -3.8852281490,
          -13.7241916520,   -12.0132825337,    -5.7025344378,    -5.0347701324,
          -11.8429453263,   -10.3608882887,    -4.9060088599,    -4.3292441889,
          -24.6981285518,   -21.5803044464,   -10.0217224305,    -8.8335740296,
          -20.8216209716,   -18.1884516718,    -8.4398899448,    -7.4373564684,
          -26.2658338232,   -22.9470218625,   -10.6387492866,    -9.3762603006,
          -22.1042889210,   -19.3066750122,    -8.9447300998,    -7.8813725083,
          -37.2397707231,   -32.5140437752,   -14.9579372792,   -13.1750641977,
          -31.0829645663,   -27.1342383953,   -12.4786111847,   -10.9894847878,
          -38.8074759945,   -33.8807611913,   -15.5749641353,   -13.7177504688,
          -32.3656325156,   -28.2524617358,   -12.9834513397,   -11.4335008277,
           -6.2697530864,    -5.4941911997,    -2.6760177087,    -2.3646979235,
           -5.5514590348,    -4.8593527871,    -2.3540418765,    -2.0780372484,
           -7.3671467764,    -6.4508933910,    -3.1079365079,    -2.7445783133,
           -6.4493265993,    -5.6421091254,    -2.7074299850,    -2.3888484763,
          -15.0489026063,   -13.1478087299,    -6.1313681028,    -5.4037410413,
          -12.7343995511,   -11.1214034936,    -5.1811467444,    -4.5645270720,
          -16.1462962963,   -14.1045109212,    -6.5632869021,    -5.7836214310,
          -13.6322671156,   -11.9041598319,    -5.5345348529,    -4.8753382999,
          -23.8280521262,   -20.8014262601,    -9.5867184969,    -8.4427841590,
          -19.9173400673,   -17.3834542001,    -8.0082516123,    -7.0510168956,
          -24.9254458162,   -21.7581284513,   -10.0186372962,    -8.8226645487,
          -20.8152076319,   -18.1662105384,    -8.3616397208,    -7.3618281235,
          -32.6072016461,   -28.4550437902,   -13.0420688910,   -11.4818272767,
          -27.1002805836,   -23.6455049066,   -10.8353564803,    -9.5375067192,
          -33.7045953361,   -29.4117459815,   -13.4739876903,   -11.8617076664,
          -27.9981481481,   -24.4282612449,   -11.1887445887,    -9.8483179471,
            0.4243484225,     0.3100365798,    -0.1852859662,    -0.1868037532,
           -0.3737560793,    -0.3713579288,    -0.4339664870,    -0.3997275789,
           -1.4046410608,    -1.2844670723,    -0.9051506317,    -0.8199377360,
           -1.8702020202,    -1.6759518260,    -1.0229466678,    -0.9177462922,
           -3.2336305441,    -2.8789707244,    -1.6250152971,    -1.4530717189,
           -3.3666479611,    -2.9805457231,    -1.6119268487,    -1.4357650054,
           -5.0626200274,    -4.4734743766,    -2.3448799626,    -2.0862057018,
           -4.8630939020,    -4.2851396203,    -2.2009070295,    -1.9537837187,
           -6.8916095107,    -6.0679780287,    -3.0647446280,    -2.7193396846,
           -6.3595398429,    -5.5897335175,    -2.7898872103,    -2.4718024319,
           -8.7205989941,    -7.6624816808,    -3.7846092935,    -3.3524736675,
           -7.8559857838,    -6.8943274147,    -3.3788673911,    -2.9898211452,
          -10.5495884774,    -9.2569853329,    -4.5044739589,    -3.9856076503,
           -9.3524317247,    -8.1989213119,    -3.9678475719,    -3.5078398584,
          -12.3785779607,   -10.8514889850,    -5.2243386243,    -4.6187416332,
          -10.8488776655,    -9.5035152091,    -4.5568277527,    -4.0258585717,
          -28.9309327846,   -25.2553673807,   -11.6191364024,   -10.2324573390,
          -24.1423057738,   -21.0710743075,    -9.6907716622,    -8.5325622424,
          -30.1502591068,   -26.3183698154,   -12.0990461793,   -10.6545466609,
          -25.1399364010,   -21.9408035723,   -10.0834251161,    -8.8779080512,
          -31.3695854291,   -27.3813722502,   -12.5789559563,   -11.0766359828,
          -26.1375670283,   -22.8105328371,   -10.4760785700,    -9.2232538600,
          -32.5889117513,   -28.4443746849,   -13.0588657332,   -11.4987253047,
          -27.1351976556,   -23.6802621019,   -10.8687320238,    -9.5685996689,
          -33.8082380735,   -29.5073771197,   -13.5387755102,   -11.9208146266,
          -28.1328282828,   -24.5499913667,   -11.2613854777,    -9.9139454777,
          -35.0275643957,   -30.5703795544,   -14.0186852872,   -12.3429039485,
          -29.1304589101,   -25.4197206314,   -11.6540389316,   -10.2592912865,
          -36.2468907179,   -31.6333819891,   -14.4985950641,   -12.7649932704,
          -30.1280895373,   -26.2894498962,   -12.0466923855,   -10.6046370954,
          -37.4662170401,   -32.6963844239,   -14.9785048411,   -13.1870825923,
          -31.1257201646,   -27.1591791610,   -12.4393458393,   -10.9499829042};

    double eps = 1.0e-6;
    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 2; ++n)
                {
                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + indexer->l2l(i, j, k) * 4;
                        int global = g + n * 2 + indexer->l2g(i, j, k) * 4;
                        EXPECT_EQ(local, system->index(
                                      g, n, indexer->l2l(i, j, k)));
                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }

    // get the boundary indexers for each face
    RCP_Bnd_Indexer lox = system->bnd_indexer(0);
    RCP_Bnd_Indexer hix = system->bnd_indexer(1);
    RCP_Bnd_Indexer loy = system->bnd_indexer(2);
    RCP_Bnd_Indexer hiy = system->bnd_indexer(3);
    RCP_Bnd_Indexer loz = system->bnd_indexer(4);
    RCP_Bnd_Indexer hiz = system->bnd_indexer(5);

    if (nodes == 1)
    {
        EXPECT_FALSE(lox.is_null());
        EXPECT_FALSE(loy.is_null());
        EXPECT_FALSE(loz.is_null());
        EXPECT_FALSE(hix.is_null());
        EXPECT_FALSE(hiy.is_null());
        EXPECT_FALSE(hiz.is_null());
    }

    if (nodes == 2)
    {
        if (node == 0)
        {
            EXPECT_FALSE(lox.is_null());
            EXPECT_FALSE(loy.is_null());
            EXPECT_FALSE(loz.is_null());
            EXPECT_TRUE(hix.is_null());
            EXPECT_FALSE(hiy.is_null());
            EXPECT_FALSE(hiz.is_null());
        }
        if (node == 1)
        {
            EXPECT_TRUE(lox.is_null());
            EXPECT_FALSE(loy.is_null());
            EXPECT_FALSE(loz.is_null());
            EXPECT_FALSE(hix.is_null());
            EXPECT_FALSE(hiy.is_null());
            EXPECT_FALSE(hiz.is_null());
        }
    }

    if (nodes == 4)
    {
        if (node == 0)
        {
            EXPECT_FALSE(lox.is_null());
            EXPECT_FALSE(loy.is_null());
            EXPECT_FALSE(loz.is_null());
            EXPECT_TRUE(hix.is_null());
            EXPECT_TRUE(hiy.is_null());
            EXPECT_FALSE(hiz.is_null());
        }
        if (node == 1)
        {
            EXPECT_TRUE(lox.is_null());
            EXPECT_FALSE(loy.is_null());
            EXPECT_FALSE(loz.is_null());
            EXPECT_FALSE(hix.is_null());
            EXPECT_TRUE(hiy.is_null());
            EXPECT_FALSE(hiz.is_null());
        }
        if (node == 2)
        {
            EXPECT_FALSE(lox.is_null());
            EXPECT_TRUE(loy.is_null());
            EXPECT_FALSE(loz.is_null());
            EXPECT_TRUE(hix.is_null());
            EXPECT_FALSE(hiy.is_null());
            EXPECT_FALSE(hiz.is_null());
        }
        if (node == 3)
        {
            EXPECT_TRUE(lox.is_null());
            EXPECT_TRUE(loy.is_null());
            EXPECT_FALSE(loz.is_null());
            EXPECT_FALSE(hix.is_null());
            EXPECT_FALSE(hiy.is_null());
            EXPECT_FALSE(hiz.is_null());
        }
    }

    if (!lox.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                for (int n = 0; n < 2; ++n)
                {
                    int lf = system->vol_unknowns()/4 + lox->local(j, k);
                    int gf = 4*4*4 + lox->l2g(j, k);

                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + lf * 4;
                        int global = g + n * 2 + gf * 4;

                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }

    if (!hix.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                for (int n = 0; n < 2; ++n)
                {
                    int lf = system->vol_unknowns()/4 + hix->local(j, k);
                    int gf = 4*4*4 + hix->l2g(j, k);

                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + lf * 4;
                        int global = g + n * 2 + gf * 4;

                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }

    if (!loy.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 2; ++n)
                {
                    int lf = system->vol_unknowns()/4 + loy->local(i, k);
                    int gf = 4*4*4 + loy->l2g(i, k);

                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + lf * 4;
                        int global = g + n * 2 + gf * 4;

                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }

    if (!hiy.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 2; ++n)
                {
                    int lf = system->vol_unknowns()/4 + hiy->local(i, k);
                    int gf = 4*4*4 + hiy->l2g(i, k);

                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + lf * 4;
                        int global = g + n * 2 + gf * 4;

                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }

    if (!loz.is_null())
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 2; ++n)
                {
                    int lf = system->vol_unknowns()/4 + loz->local(i, j);
                    int gf = 4*4*4 + loz->l2g(i, j);

                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + lf * 4;
                        int global = g + n * 2 + gf * 4;

                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }

    if (!hiz.is_null())
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 2; ++n)
                {
                    int lf = system->vol_unknowns()/4 + hiz->local(i, j);
                    int gf = 4*4*4 + hiz->l2g(i, j);

                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 2 + lf * 4;
                        int global = g + n * 2 + gf * 4;

                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP7_3Grp_Refl_RHS)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;

    using def::I; using def::J; using def::K;

    this->build(7, 3);
    this->make_data();
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    EXPECT_EQ(4, dim->num_equations());

    // make the source
    External_Source q(mesh->num_cells());
    {
        External_Source::Source_Shapes s(3, External_Source::Shape(3, 0.0));

        // source 0
        s[0][0] = 0.1;
        s[0][1] = 0.2;
        s[0][2] = 0.3;

        // source 1
        s[1][0] = 1.1;
        s[1][1] = 1.2;
        s[1][2] = 1.3;

        // source 2 is a zero source

        // source 0 -> global cells (i = 1, j,k)
        // source 1 -> global cells (i = 2, j,k)
        // all other cells have 0 source
        External_Source::ID_Field srcids(4*4*4, 2);
        for (int k = 0; k < 4; ++k)
        {
            for (int j = 0; j < 4; ++j)
            {
                srcids[indexer->g2g(1, j, k)] = 0;
                srcids[indexer->g2g(2, j, k)] = 1;
            }
        }

        // local source ids
        External_Source::ID_Field local_ids(mesh->num_cells(), 2);
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_dim(I); ++i)
                {
                    local_ids[indexer->l2l(i, j, k)] =
                             srcids[indexer->l2g(i, j, k)];
                }
            }
        }

        External_Source::Source_Field source(mesh->num_cells(), 1.0);

        q.set(local_ids, s, source);
    }

    system->build_RHS(q);

    // check q
    Teuchos::RCP<const Vector_t> rhs = system->get_RHS();
    EXPECT_EQ(mesh->num_cells() * 4 * 3, VectorTraits::local_length(rhs));
    Teuchos::ArrayView<const double> rhs_data = VectorTraits::get_data(rhs,0);

    double eps = 1.0e-6;
    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int g_i  = i + indexer->offset(I);
                int cell = indexer->l2l(i, j, k);
                double f = 1.0;

                if (g_i == 1)
                {
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 0 + cell * 12], eps);

                    f = -2.0/3.0;
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 1 + cell * 12], eps);

                    f = 8.0/15.0;
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 2 + cell * 12], eps);

                    f = -16.0/35.0;
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 3 + cell * 12], eps);
                }
                else if (g_i == 2)
                {
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 0 + cell * 12], eps);

                    f = -2.0/3.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 1 + cell * 12], eps);

                    f = 8.0/15.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 2 + cell * 12], eps);

                    f = -16.0/35.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 3 + cell * 12], eps);
                }
                else
                {
                    for (int n = 0; n < 4; ++n)
                    {
                        EXPECT_SOFTEQ(0.0, rhs_data[0 + 3 * n + cell * 12], eps);
                        EXPECT_SOFTEQ(0.0, rhs_data[1 + 3 * n + cell * 12], eps);
                        EXPECT_SOFTEQ(0.0, rhs_data[2 + 3 * n + cell * 12], eps);
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP7_3Grp_Isotropic_RHS)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;

    using def::I; using def::J; using def::K;
    this->build(7, 3);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    db->set("boundary", string("isotropic"));
    {
        Partitioner::ParameterList bdb("boundary");
        Array_Dbl lx(3, 1.1);
        Array_Dbl hx(3, 1.2);
        Array_Dbl ly(3, 1.3);
        Array_Dbl hy(3, 1.4);
        Array_Dbl lz(3, 1.5);
        Array_Dbl hz(3, 1.6);

        bdb.set("minus_x_phi", lx);
        bdb.set("plus_x_phi", hx);
        bdb.set("minus_y_phi", ly);
        bdb.set("plus_y_phi", hy);
        bdb.set("minus_z_phi", lz);
        bdb.set("plus_z_phi", hz);

        db->set("boundary_db", bdb);
    }

    this->make_data();
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    int num_groups = this->num_groups;
    EXPECT_EQ(4, dim->num_equations());

    // make the source
    External_Source q(mesh->num_cells());
    {
        External_Source::Source_Shapes s(3, External_Source::Shape(3, 0.0));

        // source 0
        s[0][0] = 0.1;
        s[0][1] = 0.2;
        s[0][2] = 0.3;

        // source 1
        s[1][0] = 1.1;
        s[1][1] = 1.2;
        s[1][2] = 1.3;

        // source 2 is a zero source

        // source 0 -> global cells (i = 1, j,k)
        // source 1 -> global cells (i = 2, j,k)
        // all other cells have 0 source
        External_Source::ID_Field srcids(4*4*4, 2);
        for (int k = 0; k < 4; ++k)
        {
            for (int j = 0; j < 4; ++j)
            {
                srcids[indexer->g2g(1, j, k)] = 0;
                srcids[indexer->g2g(2, j, k)] = 1;
            }
        }

        // local source ids
        External_Source::ID_Field local_ids(mesh->num_cells(), 2);
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_dim(I); ++i)
                {
                    local_ids[indexer->l2l(i, j, k)] =
                             srcids[indexer->l2g(i, j, k)];
                }
            }
        }

        External_Source::Source_Field source(mesh->num_cells(), 1.0);

        q.set(local_ids, s, source);
    }

    system->build_RHS(q);

    // check q
    Teuchos::RCP<const Vector_t> rhs = system->get_RHS();
    EXPECT_EQ(mesh->num_cells() * 4 * 3, system->vol_unknowns());
    EXPECT_EQ((system->vol_unknowns() + system->bnd_unknowns()),
              VectorTraits::local_length(rhs));
    Teuchos::ArrayView<const double> rhs_data = VectorTraits::get_data(rhs,0);

    double eps = 1.0e-6;
    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int g_i  = i + indexer->offset(I);
                int cell = indexer->l2l(i, j, k);
                double f = 1.0;

                if (g_i == 1)
                {
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 0 + cell * 12], eps);

                    f = -2.0/3.0;
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 1 + cell * 12], eps);

                    f = 8.0/15.0;
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 2 + cell * 12], eps);

                    f = -16.0/35.0;
                    EXPECT_SOFTEQ(f*0.1, rhs_data[0 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.2, rhs_data[1 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*0.3, rhs_data[2 + 3 * 3 + cell * 12], eps);
                }
                else if (g_i == 2)
                {
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 0 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 0 + cell * 12], eps);

                    f = -2.0/3.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 1 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 1 + cell * 12], eps);

                    f = 8.0/15.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 2 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 2 + cell * 12], eps);

                    f = -16.0/35.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[0 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.2, rhs_data[1 + 3 * 3 + cell * 12], eps);
                    EXPECT_SOFTEQ(f*1.3, rhs_data[2 + 3 * 3 + cell * 12], eps);
                }
                else
                {
                    for (int n = 0; n < 4; ++n)
                    {
                        EXPECT_SOFTEQ(0.0, rhs_data[0 + 3 * n + cell * 12], eps);
                        EXPECT_SOFTEQ(0.0, rhs_data[1 + 3 * n + cell * 12], eps);
                        EXPECT_SOFTEQ(0.0, rhs_data[2 + 3 * n + cell * 12], eps);
                    }
                }
            }
        }
    }

    // check faces

    // get the boundary indexers for each face
    RCP_Bnd_Indexer lox = system->bnd_indexer(0);
    RCP_Bnd_Indexer hix = system->bnd_indexer(1);
    RCP_Bnd_Indexer loy = system->bnd_indexer(2);
    RCP_Bnd_Indexer hiy = system->bnd_indexer(3);
    RCP_Bnd_Indexer loz = system->bnd_indexer(4);
    RCP_Bnd_Indexer hiz = system->bnd_indexer(5);

    double f = 0.0;

    if (!lox.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                int lf = system->vol_unknowns()/12 + lox->local(j, k);

                for (int g = 0; g < num_groups; ++g)
                {
                    f = 1.0/2.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[g + 3 * 0 + lf * 12], eps);

                    f = -1.0/8.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[g + 3 * 1 + lf * 12], eps);

                    f = 1.0/16.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[g + 3 * 2 + lf * 12], eps);

                    f = -5.0/128.0;
                    EXPECT_SOFTEQ(f*1.1, rhs_data[g + 3 * 3 + lf * 12], eps);
                }
            }
        }
    }

    if (!hix.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                int lf = system->vol_unknowns()/12 + hix->local(j, k);

                for (int g = 0; g < num_groups; ++g)
                {
                    f = 1.0/2.0;
                    EXPECT_SOFTEQ(f*1.2, rhs_data[g + 3 * 0 + lf * 12], eps);

                    f = -1.0/8.0;
                    EXPECT_SOFTEQ(f*1.2, rhs_data[g + 3 * 1 + lf * 12], eps);

                    f = 1.0/16.0;
                    EXPECT_SOFTEQ(f*1.2, rhs_data[g + 3 * 2 + lf * 12], eps);

                    f = -5.0/128.0;
                    EXPECT_SOFTEQ(f*1.2, rhs_data[g + 3 * 3 + lf * 12], eps);
                }
            }
        }
    }

    if (!loy.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int lf = system->vol_unknowns()/12 + loy->local(i, k);

                for (int g = 0; g < num_groups; ++g)
                {
                    f = 1.0/2.0;
                    EXPECT_SOFTEQ(f*1.3, rhs_data[g + 3 * 0 + lf * 12], eps);

                    f = -1.0/8.0;
                    EXPECT_SOFTEQ(f*1.3, rhs_data[g + 3 * 1 + lf * 12], eps);

                    f = 1.0/16.0;
                    EXPECT_SOFTEQ(f*1.3, rhs_data[g + 3 * 2 + lf * 12], eps);

                    f = -5.0/128.0;
                    EXPECT_SOFTEQ(f*1.3, rhs_data[g + 3 * 3 + lf * 12], eps);
                }
            }
        }
    }

    if (!hiy.is_null())
    {
        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int lf = system->vol_unknowns()/12 + hiy->local(i, k);

                for (int g = 0; g < num_groups; ++g)
                {
                    f = 1.0/2.0;
                    EXPECT_SOFTEQ(f*1.4, rhs_data[g + 3 * 0 + lf * 12], eps);

                    f = -1.0/8.0;
                    EXPECT_SOFTEQ(f*1.4, rhs_data[g + 3 * 1 + lf * 12], eps);

                    f = 1.0/16.0;
                    EXPECT_SOFTEQ(f*1.4, rhs_data[g + 3 * 2 + lf * 12], eps);

                    f = -5.0/128.0;
                    EXPECT_SOFTEQ(f*1.4, rhs_data[g + 3 * 3 + lf * 12], eps);
                }
            }
        }
    }

    for (int j = 0; j < mesh->num_cells_dim(J); ++j)
    {
        for (int i = 0; i < mesh->num_cells_dim(I); ++i)
        {
            int lf = system->vol_unknowns()/12 + loz->local(i, j);

            for (int g = 0; g < num_groups; ++g)
            {
                f = 1.0/2.0;
                EXPECT_SOFTEQ(f*1.5, rhs_data[g + 3 * 0 + lf * 12], eps);

                f = -1.0/8.0;
                EXPECT_SOFTEQ(f*1.5, rhs_data[g + 3 * 1 + lf * 12], eps);

                f = 1.0/16.0;
                EXPECT_SOFTEQ(f*1.5, rhs_data[g + 3 * 2 + lf * 12], eps);

                f = -5.0/128.0;
                EXPECT_SOFTEQ(f*1.5, rhs_data[g + 3 * 3 + lf * 12], eps);
            }
        }
    }

    for (int j = 0; j < mesh->num_cells_dim(J); ++j)
    {
        for (int i = 0; i < mesh->num_cells_dim(I); ++i)
        {
            int lf = system->vol_unknowns()/12 + hiz->local(i, j);

            for (int g = 0; g < num_groups; ++g)
            {
                f = 1.0/2.0;
                EXPECT_SOFTEQ(f*1.6, rhs_data[g + 3 * 0 + lf * 12], eps);

                f = -1.0/8.0;
                EXPECT_SOFTEQ(f*1.6, rhs_data[g + 3 * 1 + lf * 12], eps);

                f = 1.0/16.0;
                EXPECT_SOFTEQ(f*1.6, rhs_data[g + 3 * 2 + lf * 12], eps);

                f = -5.0/128.0;
                EXPECT_SOFTEQ(f*1.6, rhs_data[g + 3 * 3 + lf * 12], eps);
            }
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP7_3Grp_Null_Fission_Matrix)
{
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;
    typedef typename TypeParam::MV                  MV;
    typedef typename TypeParam::OP                  OP;
    typedef Anasazi::OperatorTraits<double,MV,OP>   OPT;

    using def::I; using def::J; using def::K;
    this->build(7, 3);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    // 2 materials-no fission
    vector<int>    ids(2, 0);
    vector<double> f(2, 0.0);
    vector<int>    matids(mesh->num_cells(), 0);
    ids[0] = 9;  f[0] = 1.0;
    ids[1] = 11; f[1] = 1.0;

    vector<int> gids(4*4*4, 0);
    for (int k = 0; k < 4; ++k)
    {
        for (int j = 0; j < 4; ++j)
        {
            gids[indexer->g2g(0, j, k)] = 9;
            gids[indexer->g2g(1, j, k)] = 11;
            gids[indexer->g2g(2, j, k)] = 9;
            gids[indexer->g2g(3, j, k)] = 11;
        }
    }

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                matids[indexer->l2l(i, j, k)] = gids[indexer->l2g(i, j, k)];
            }
        }
    }

    this->make_data(ids, f, matids);
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    int num_groups = this->num_groups;
    EXPECT_EQ(4, dim->num_equations());

    // make the matrix
    system->build_fission_matrix();
    Teuchos::RCP<const Matrix_t> B =
        Teuchos::rcp_dynamic_cast<const Matrix_t>(system->get_fission_matrix());

    EXPECT_EQ(12 * mesh->num_cells(), MatrixTraits::local_rows(B));
    EXPECT_EQ(12 * data->num_cells(), MatrixTraits::global_rows(B));
    EXPECT_EQ(12 * data->num_cells(), MatrixTraits::global_columns(B));

    // it should be all zeros

    // make a vector
    Teuchos::RCP<Vector_t> x = VectorTraits::build_vector(system->get_Map());
    Teuchos::RCP<Vector_t> y = VectorTraits::build_vector(system->get_Map());
    EXPECT_EQ(mesh->num_cells() * 12, VectorTraits::local_length(x));
    EXPECT_EQ(mesh->num_cells() * 12, VectorTraits::local_length(y));
    Teuchos::ArrayView<double> x_data = VectorTraits::get_data_nonconst(x,0);
    Teuchos::ArrayView<double> y_data = VectorTraits::get_data_nonconst(y,0);
    VectorTraits::put_scalar(y,-1.0);

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int global = indexer->l2g(i, j, k);
                int local  = indexer->l2l(i, j, k);

                for (int n = 0; n < 2; ++n)
                {
                    double gv = static_cast<double>(global + n);
                    x_data[0 + n * num_groups + local * num_groups * 2] = gv + 0.1;
                    x_data[1 + n * num_groups + local * num_groups * 2] = gv + 0.2;
                }
            }
        }
    }

    OPT::Apply(*B,*x,*y);

    for (int n = 0; n < y_data.size(); ++n)
    {
        EXPECT_EQ(0.0, y_data[n]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(MatrixTest, SP7_3Grp_Fission_Matrix)
{
    typedef typename TestFixture::Linear_System     Linear_System;
    typedef typename TestFixture::RCP_Linear_System RCP_Linear_System;
    typedef typename TestFixture::Matrix_t          Matrix_t;
    typedef typename TestFixture::Vector_t          Vector_t;
    typedef profugus::MatrixTraits<TypeParam>       MatrixTraits;
    typedef profugus::VectorTraits<TypeParam>       VectorTraits;
    typedef typename TypeParam::MV                  MV;
    typedef typename TypeParam::OP                  OP;
    typedef Anasazi::OperatorTraits<double,MV,OP>   OPT;

    using def::I; using def::J; using def::K;
    this->build(7, 3);
    RCP_ParameterList db = this->db;
    RCP_Mesh mesh = this->mesh;
    RCP_Indexer indexer = this->indexer;
    RCP_Global_Data data = this->data;

    // 2 materials (mat 9 has fission)
    vector<int>    ids(2, 0);
    vector<double> f(2, 0.0);
    vector<int>    matids(mesh->num_cells(), 0);
    ids[0] = 9;  f[0] = 1.0;
    ids[1] = 11; f[1] = 1.0;

    vector<int> gids(4*4*4, 0);
    for (int k = 0; k < 4; ++k)
    {
        for (int j = 0; j < 4; ++j)
        {
            gids[indexer->g2g(0, j, k)] = 9;
            gids[indexer->g2g(1, j, k)] = 11;
            gids[indexer->g2g(2, j, k)] = 9;
            gids[indexer->g2g(3, j, k)] = 11;
        }
    }

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                matids[indexer->l2l(i, j, k)] = gids[indexer->l2g(i, j, k)];
            }
        }
    }

    this->make_data(ids, f, matids);
    RCP_Dimensions dim = this->dim;
    RCP_Linear_System system = this->system;
    int num_groups = this->num_groups;
    EXPECT_EQ(4, dim->num_equations());

    // add fission to mat (we don't need to make the scattering matrices as
    // they don't intrude on the Fission matrix)
    RCP_Mat_DB matf = Teuchos::rcp(new Mat_DB_t);
    {
        const XS &old = this->mat->xs();
        RCP_XS xs     = Teuchos::rcp(new XS);

        xs->set(old.pn_order(), old.num_groups());

        XS::OneDArray tot(3), nusigf(3), chi(3);
        XS::TwoDArray P0(3, 3), P1(3, 3), P2(3, 3), P3(3, 3);

        // material 9
        {
            const XS::Vector &sig = old.vector(9, XS::TOTAL);
            for (int g = 0; g < 3; ++g)
            {
                tot[g] = sig[g];
            }
            nusigf[0] = 1.3;
            nusigf[1] = 4.2;
            nusigf[2] = 0.0;

            chi[0] = 0.2;
            chi[1] = 0.8;
            chi[2] = 0.0;

            xs->add(9, XS::TOTAL, tot);
            xs->add(9, XS::NU_SIG_F, nusigf);
            xs->add(9, XS::CHI, chi);
        }

        // material 11
        {
            const XS::Vector &sig = old.vector(11, XS::TOTAL);
            for (int g = 0; g < 3; ++g)
            {
                tot[g] = sig[g];
            }

            xs->add(11, XS::TOTAL, tot);
        }

        xs->complete();
        matf->set(xs, mesh->num_cells());

        for (int k = 0; k < mesh->num_cells_dim(K); ++k)
        {
            for (int j = 0; j < mesh->num_cells_dim(J); ++j)
            {
                for (int i = 0; i < mesh->num_cells_dim(I); ++i)
                {
                    matf->matid(indexer->l2l(i ,j, k)) =
                        gids[indexer->l2g(i, j, k)];
                }
            }
        }

        system = Teuchos::rcp(new Linear_System(
                                  db, dim, matf, mesh, indexer, data));
    }

    // make the matrix
    system->build_fission_matrix();
    Teuchos::RCP<const Matrix_t> B =
        Teuchos::rcp_dynamic_cast<const Matrix_t>(system->get_fission_matrix());

    EXPECT_EQ(12 * mesh->num_cells(), MatrixTraits::local_rows(B));
    EXPECT_EQ(12 * data->num_cells(), MatrixTraits::global_rows(B));
    EXPECT_EQ(12 * data->num_cells(), MatrixTraits::global_columns(B));

    // extract some blocks
    if (node == 0)
    {
        // 2nd equation, 2nd cell
        int row = 28;

        Teuchos::ArrayView<const int>    indices;
        Teuchos::ArrayView<const double> values;
        MatrixTraits::get_local_row_view(B,row,indices,values);
        EXPECT_EQ(8, indices.size());
        EXPECT_EQ(8, values.size());

        if (nodes == 1)
        {
            EXPECT_EQ(24, MatrixTraits::global_col_id(B,indices[0]));
            EXPECT_EQ(25, MatrixTraits::global_col_id(B,indices[1]));
            EXPECT_EQ(27, MatrixTraits::global_col_id(B,indices[2]));
            EXPECT_EQ(28, MatrixTraits::global_col_id(B,indices[3]));
            EXPECT_EQ(30, MatrixTraits::global_col_id(B,indices[4]));
            EXPECT_EQ(31, MatrixTraits::global_col_id(B,indices[5]));
            EXPECT_EQ(33, MatrixTraits::global_col_id(B,indices[6]));
            EXPECT_EQ(34, MatrixTraits::global_col_id(B,indices[7]));

            EXPECT_SOFTEQ(-0.69333333, values[0], 1.0e-6);
            EXPECT_SOFTEQ(-2.24000000, values[1], 1.0e-6);
            EXPECT_SOFTEQ( 0.46222222, values[2], 1.0e-6);
            EXPECT_SOFTEQ( 1.49333333, values[3], 1.0e-6);
            EXPECT_SOFTEQ(-0.36977778, values[4], 1.0e-6);
            EXPECT_SOFTEQ(-1.19466667, values[5], 1.0e-6);
            EXPECT_SOFTEQ( 0.31695238, values[6], 1.0e-6);
            EXPECT_SOFTEQ( 1.02400000, values[7], 1.0e-6);
        }
    }

    // make a vector
    Teuchos::RCP<Vector_t> x = VectorTraits::build_vector(system->get_Map());
    Teuchos::RCP<Vector_t> y = VectorTraits::build_vector(system->get_Map());
    EXPECT_EQ(mesh->num_cells() * 12, VectorTraits::local_length(x));
    EXPECT_EQ(mesh->num_cells() * 12, VectorTraits::local_length(y));
    Teuchos::ArrayView<double> x_data = VectorTraits::get_data_nonconst(x,0);
    Teuchos::ArrayView<double> y_data = VectorTraits::get_data_nonconst(y,0);
    VectorTraits::put_scalar(y,-1.0);

    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                int global = indexer->l2g(i, j, k);
                int local  = indexer->l2l(i, j, k);

                for (int n = 0; n < 4; ++n)
                {
                    double gv = static_cast<double>(global + n);
                    x_data[0 + n * num_groups + local * num_groups * 4] = gv + 0.1;
                    x_data[1 + n * num_groups + local * num_groups * 4] = gv + 0.2;
                    x_data[2 + n * num_groups + local * num_groups * 4] = gv + 0.3;
                }
            }
        }
    }

    OPT::Apply(*B,*x,*y);

    double v[768] = {
        -0.989123810,   -3.956495238,    0.000000000,    0.659415873,
        2.637663492,    0.000000000,   -0.527532698,   -2.110130794,
        0.000000000,    0.452170884,    1.808683537,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        -0.088171429,   -0.352685714,    0.000000000,    0.058780952,
        0.235123810,    0.000000000,   -0.047024762,   -0.188099048,
        0.000000000,    0.040306939,    0.161227755,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.812780952,    3.251123810,    0.000000000,   -0.541853968,
        -2.167415873,    0.000000000,    0.433483175,    1.733932698,
        0.000000000,   -0.371557007,   -1.486228027,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        1.713733333,    6.854933333,    0.000000000,   -1.142488889,
        -4.569955556,    0.000000000,    0.913991111,    3.655964444,
        0.000000000,   -0.783420952,   -3.133683810,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        2.614685714,   10.458742857,    0.000000000,   -1.743123810,
        -6.972495238,    0.000000000,    1.394499048,    5.577996190,
        0.000000000,   -1.195284898,   -4.781139592,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        3.515638095,   14.062552381,    0.000000000,   -2.343758730,
        -9.375034921,    0.000000000,    1.875006984,    7.500027937,
        0.000000000,   -1.607148844,   -6.428595374,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        4.416590476,   17.666361905,    0.000000000,   -2.944393651,
        -11.777574603,    0.000000000,    2.355514921,    9.422059683,
        0.000000000,   -2.019012789,   -8.076051156,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        5.317542857,   21.270171429,    0.000000000,   -3.545028571,
        -14.180114286,    0.000000000,    2.836022857,   11.344091429,
        0.000000000,   -2.430876735,   -9.723506939,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        6.218495238,   24.873980952,    0.000000000,   -4.145663492,
        -16.582653968,    0.000000000,    3.316530794,   13.266123175,
        0.000000000,   -2.842740680,  -11.370962721,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        7.119447619,   28.477790476,    0.000000000,   -4.746298413,
        -18.985193651,    0.000000000,    3.797038730,   15.188154921,
        0.000000000,   -3.254604626,  -13.018418503,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        8.020400000,   32.081600000,    0.000000000,   -5.346933333,
        -21.387733333,    0.000000000,    4.277546667,   17.110186667,
        0.000000000,   -3.666468571,  -14.665874286,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        8.921352381,   35.685409524,    0.000000000,   -5.947568254,
        -23.790273016,    0.000000000,    4.758054603,   19.032218413,
        0.000000000,   -4.078332517,  -16.313330068,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        9.822304762,   39.289219048,    0.000000000,   -6.548203175,
        -26.192812698,    0.000000000,    5.238562540,   20.954250159,
        0.000000000,   -4.490196463,  -17.960785850,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        10.723257143,   42.893028571,    0.000000000,   -7.148838095,
        -28.595352381,    0.000000000,    5.719070476,   22.876281905,
        0.000000000,   -4.902060408,  -19.608241633,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        11.624209524,   46.496838095,    0.000000000,   -7.749473016,
        -30.997892063,    0.000000000,    6.199578413,   24.798313651,
        0.000000000,   -5.313924354,  -21.255697415,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        12.525161905,   50.100647619,    0.000000000,   -8.350107937,
        -33.400431746,    0.000000000,    6.680086349,   26.720345397,
        0.000000000,   -5.725788299,  -22.903153197,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        13.426114286,   53.704457143,    0.000000000,   -8.950742857,
        -35.802971429,    0.000000000,    7.160594286,   28.642377143,
        0.000000000,   -6.137652245,  -24.550608980,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        14.327066667,   57.308266667,    0.000000000,   -9.551377778,
        -38.205511111,    0.000000000,    7.641102222,   30.564408889,
        0.000000000,   -6.549516190,  -26.198064762,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        15.228019048,   60.912076190,    0.000000000,  -10.152012698,
        -40.608050794,    0.000000000,    8.121610159,   32.486440635,
        0.000000000,   -6.961380136,  -27.845520544,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        16.128971429,   64.515885714,    0.000000000,  -10.752647619,
        -43.010590476,    0.000000000,    8.602118095,   34.408472381,
        0.000000000,   -7.373244082,  -29.492976327,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        17.029923810,   68.119695238,    0.000000000,  -11.353282540,
        -45.413130159,    0.000000000,    9.082626032,   36.330504127,
        0.000000000,   -7.785108027,  -31.140432109,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        17.930876190,   71.723504762,    0.000000000,  -11.953917460,
        -47.815669841,    0.000000000,    9.563133968,   38.252535873,
        0.000000000,   -8.196971973,  -32.787887891,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        18.831828571,   75.327314286,    0.000000000,  -12.554552381,
        -50.218209524,    0.000000000,   10.043641905,   40.174567619,
        0.000000000,   -8.608835918,  -34.435343673,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        19.732780952,   78.931123810,    0.000000000,  -13.155187302,
        -52.620749206,    0.000000000,   10.524149841,   42.096599365,
        0.000000000,   -9.020699864,  -36.082799456,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        20.633733333,   82.534933333,    0.000000000,  -13.755822222,
        -55.023288889,    0.000000000,   11.004657778,   44.018631111,
        0.000000000,   -9.432563810,  -37.730255238,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        21.534685714,   86.138742857,    0.000000000,  -14.356457143,
        -57.425828571,    0.000000000,   11.485165714,   45.940662857,
        0.000000000,   -9.844427755,  -39.377711020,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        22.435638095,   89.742552381,    0.000000000,  -14.957092063,
        -59.828368254,    0.000000000,   11.965673651,   47.862694603,
        0.000000000,  -10.256291701,  -41.025166803,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        23.336590476,   93.346361905,    0.000000000,  -15.557726984,
        -62.230907937,    0.000000000,   12.446181587,   49.784726349,
        0.000000000,  -10.668155646,  -42.672622585,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        24.237542857,   96.950171429,    0.000000000,  -16.158361905,
        -64.633447619,    0.000000000,   12.926689524,   51.706758095,
        0.000000000,  -11.080019592,  -44.320078367,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        25.138495238,  100.553980952,    0.000000000,  -16.758996825,
        -67.035987302,    0.000000000,   13.407197460,   53.628789841,
        0.000000000,  -11.491883537,  -45.967534150,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        26.039447619,  104.157790476,    0.000000000,  -17.359631746,
        -69.438526984,    0.000000000,   13.887705397,   55.550821587,
        0.000000000,  -11.903747483,  -47.614989932,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        26.940400000,  107.761600000,    0.000000000,  -17.960266667,
        -71.841066667,    0.000000000,   14.368213333,   57.472853333,
        0.000000000,  -12.315611429,  -49.262445714,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
        0.000000000,    0.000000000,    0.000000000,    0.000000000,
    };

    double eps = 1.0e-5;
    for (int k = 0; k < mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < mesh->num_cells_dim(I); ++i)
            {
                for (int n = 0; n < 4; ++n)
                {
                    for (int g = 0; g < num_groups; ++g)
                    {
                        int local  = g + n * 3 + indexer->l2l(i, j, k) * 12;
                        int global = g + n * 3 + indexer->l2g(i, j, k) * 12;
                        EXPECT_EQ(local, system->index(
                                      g, n, indexer->l2l(i, j, k)));
                        EXPECT_SOFTEQ(v[global], y_data[local], eps);
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
//                 end of tstLinear_System_FV.cc
//---------------------------------------------------------------------------//
