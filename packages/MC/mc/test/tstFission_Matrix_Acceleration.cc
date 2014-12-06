//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFission_Matrix_Acceleration.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 12 14:59:19 2014
 * \brief  Test for Fission_Matrix_Acceleration
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>
#include <cmath>

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_Export.h"
#include "Teuchos_RCP.hpp"

#include <SPn/config.h>

#include "rng/RNG_Control.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "spn/Moment_Coefficients.hh"
#include "spn/MatrixTraits.hh"
#include "spn/VectorTraits.hh"
#include "../Fission_Matrix_Acceleration.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template<class T>
class FM_AccelerationTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Fission_Matrix_Acceleration         Acceleration;
    typedef Acceleration::RCP_ParameterList               RCP_ParameterList;
    typedef std::shared_ptr<Acceleration>                 SP_Acceleration;
    typedef Acceleration::Problem_Builder_t               SPN_Builder;
    typedef profugus::Fission_Matrix_Acceleration_Impl<T> Implementation;
    typedef std::shared_ptr<Implementation>               SP_Implementation;
    typedef typename Implementation::Linear_System_t      Linear_System_t;
    typedef typename Linear_System_t::RCP_Vector          RCP_Vector;
    typedef std::vector<double>                           Vec_Dbl;
    typedef Acceleration::Fission_Site                    Fission_Site;
    typedef Acceleration::Fission_Site_Container          Fission_Site_Container;
    typedef profugus::VectorTraits<T>                     Traits;
    typedef profugus::MatrixTraits<T>                     MTraits;

  protected:
    void SetUp()
    {
        // make the acceleration
        implementation =
            std::make_shared<profugus::Fission_Matrix_Acceleration_Impl<T>>();
        acceleration = implementation;

        mc_db = Teuchos::rcp(new Acceleration::ParameterList());
        mc_db->set("fission_matrix_db",
                   Acceleration::ParameterList("fission_matrix"));

        seed = 12134;
    }

    void get_adjoint(Vec_Dbl &phi)
    {
        auto u = implementation->adjoint();

        const Linear_System_t &system = implementation->spn_system();

        N  = system.get_dims()->num_equations();
        Ng = acceleration->mat().xs().num_groups();
        Nc = acceleration->mesh().num_cells();

        // SPN moments
        double u_m[4] = {0.0, 0.0, 0.0, 0.0};

        // resize phi
        phi.resize(Nc * Ng);

        // loop over groups
        for (int g = 0; g < Ng; ++g)
        {
            // loop over cells on this domain
            for (int cell = 0; cell < Nc; ++cell)
            {
                // loop over moments
                for (int n = 0; n < N; ++n)
                {
                     EXPECT_GT(u.size(), system.index(g, n, cell));

                    // get the moment from the solution vector
                    u_m[n] = u[system.index(g, n, cell)];
                }

                // assign the scalar flux (SPN Phi_0 moment)
                phi[cell + Nc * g] =  profugus::Moment_Coefficients::u_to_phi(
                    u_m[0], u_m[1], u_m[2], u_m[3]);
            }
        }
    }

    void build_fs(bool random = false)
    {
        profugus::RNG_Control rcon(seed);
        auto rng = rcon.rng();

        fs.clear();

        if (!random)
        {
            // we need 13 total sites to match the reference from the
            // Acceleration Tests notebook
            fs.resize(13);

            // 2 in cell 5
            // 4 in cell 6
            // 5 in cell 9
            // 2 in cell 10

            double z = 0.125;
            double a[] = {1.0, 2.0};
            double b[] = {2.0, 3.0};

            double x[][2] = {{1., 2.}, {1., 2.},
                             {2., 3.}, {2., 3.}, {2., 3.}, {2., 3.},
                             {1., 2.}, {1., 2.}, {1., 2.}, {1., 2.}, {1., 2.},
                             {2., 3.}, {2., 3.}};
            double y[][2] = {{1., 2.}, {1., 2.},
                             {1., 2.}, {1., 2.}, {1., 2.}, {1., 2.},
                             {2., 3.}, {2., 3.}, {2., 3.}, {2., 3.}, {2., 3.},
                             {2., 3.}, {2., 3.}};

            for (int n = 0; n < 13; ++n)
            {
                fs[n].r[0] = rng.ran() * (x[n][1] - x[n][0]) + x[n][0];
                fs[n].r[1] = rng.ran() * (y[n][1] - y[n][0]) + y[n][0];
                fs[n].r[2] = z;
            }
        }
        else
        {
            // randomly sample 20 sites in the fuel regions
            fs.resize(20);

            double z = 0.125;

            for (int n = 0; n < 20; ++n)
            {
                // sample the cell
                fs[n].r[0] = 2.0 * rng.ran() + 1.0;
                fs[n].r[1] = 2.0 * rng.ran() + 1.0;
                fs[n].r[2] = z;
            }
        }
    }

    void print_matrices()
    {
        auto A = implementation->spn_system().get_Matrix();
        auto B = Teuchos::rcp_dynamic_cast<typename T::MATRIX>(
            implementation->spn_system().get_fission_matrix());
        MTraits::write_matrix_file(A, "A.mtx");
        MTraits::write_matrix_file(B, "B.mtx");
    }

  protected:
    // >>> DATA

    RCP_ParameterList mc_db;

    SP_Implementation implementation;
    SP_Acceleration   acceleration;
    SPN_Builder       builder;

    Fission_Site_Container fs;

    int N, Ng, Nc;

    int seed;
};

//---------------------------------------------------------------------------//
// SETUP TEMPLATED TESTS
//---------------------------------------------------------------------------//
// Currently only Epetra is supported (due to adjoints)

typedef testing::Types<profugus::EpetraTypes> Test_Types;

TYPED_TEST_CASE(FM_AccelerationTest, Test_Types);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TYPED_TEST(FM_AccelerationTest, initialization)
{
    using std::string;

    // setup the spn problem
    this->builder.setup("inf_med_k.xml");

    // build the spn problem
    this->acceleration->build_problem(this->builder);

    // initialize the spn problem
    this->acceleration->initialize(this->mc_db);

    const auto &fmdb = this->mc_db->sublist("fission_matrix_db");
    EXPECT_EQ("stratimikos", fmdb.template get<string>("solver_type"));
    EXPECT_EQ("ml", fmdb.template get<string>("Preconditioner"));
    EXPECT_EQ("Belos", fmdb.sublist(
                  "Stratimikos").template get<string>("Linear Solver Type"));

    // check the reference
    typename TestFixture::Vec_Dbl adjoint;
    this->get_adjoint(adjoint);

    typename TestFixture::Vec_Dbl inf_med_vec(this->Ng, 0.0);
    double norm = 0.0;
    for (int g = 0; g < this->Ng; ++g)
    {
        for (int cell = 0; cell < this->Nc; ++cell)
        {
            inf_med_vec[g] += adjoint[cell + this->Nc * g];
        }
        inf_med_vec[g] /= static_cast<double>(this->Nc);
        norm           += inf_med_vec[g] * inf_med_vec[g];
    }
    norm = 1.0 / std::sqrt(norm);

    // normalize the vector
    for (auto &x : inf_med_vec)
    {
        x *= norm;
    }

    double ref[] = {0.33611574,  0.3155125 ,  0.31716512,  0.40852777,
                    0.41532128, 0.4177149 ,  0.41594701};

    EXPECT_VEC_SOFTEQ(ref, inf_med_vec, 1.0e-5);

    EXPECT_SOFTEQ(1.1889633277246718, this->implementation->keff(), 1.0e-6);
}

//---------------------------------------------------------------------------//

TYPED_TEST(FM_AccelerationTest, setup)
{
    using std::string;
    typedef typename TestFixture::Traits Traits;

    // setup the spn problem
    this->builder.setup("mesh4x4.xml");

    // build the spn problem
    this->acceleration->build_problem(this->builder);
    EXPECT_EQ(16, this->acceleration->mesh().num_cells());

    // initialize the spn problem
    this->acceleration->initialize(this->mc_db);

    // build fission source
    this->build_fs();

    // starting cycle
    this->acceleration->start_cycle(1.0, this->fs);

    // check the vector Bpsi
    auto Bpsi_l = this->implementation->get_Bpsi_l();

    EXPECT_EQ(16*2*3, Traits::local_length(Bpsi_l));

    // make the reference
    const auto &system = this->implementation->spn_system();

    // get the fission matrix
    auto B = system.get_fission_matrix();

    // see "Acceleration Tests" notebook for u vector that is equivalent to
    // the fission source in this test
    auto y = Traits::build_vector(Teuchos::rcpFromRef(B->OperatorRangeMap()));
    auto u = Traits::build_vector(Teuchos::rcpFromRef(B->OperatorDomainMap()));
    EXPECT_EQ(16*2*3, Traits::local_length(y));
    EXPECT_EQ(16*2*3, Traits::local_length(u));
    {
        double ref_u[] = {
            3.34382, 9.64797, 4.00988, 3.17858, 2.81054, 2.48293, 3.34382,
            9.64797, 4.00988, 3.17858, 2.81054, 2.48293, 3.34382, 9.64797,
            4.00988, 3.17858, 2.81054, 2.48293, 3.34382, 9.64797, 4.00988,
            3.17858, 2.81054, 2.48293, 3.34382, 9.64797, 4.00988, 3.17858,
            2.81054, 2.48293, 3.68000, 13.54000, 16.40000, 4.02000, 19.56000,
            24.30000, 8.78000, 6.49000, 17.80000, 8.67000, 9.36000, 22.20000,
            3.34382, 9.64797, 4.00988, 3.17858, 2.81054, 2.48293, 3.34382,
            9.64797, 4.00988, 3.17858, 2.81054, 2.48293, 3.12000, 17.42000,
            16.80000, 1.68000, 21.63000, 10.20000, 3.68000, 13.54000, 16.40000,
            4.02000, 19.56000, 24.30000, 3.34382, 9.64797, 4.00988, 3.17858,
            2.81054, 2.48293, 3.34382, 9.64797, 4.00988, 3.17858, 2.81054,
            2.48293, 3.34382, 9.64797, 4.00988, 3.17858, 2.81054, 2.48293,
            3.34382, 9.64797, 4.00988, 3.17858, 2.81054, 2.48293, 3.34382,
            9.64797, 4.00988, 3.17858, 2.81054, 2.48293};

        auto v = Traits::get_data_nonconst(u);
        std::copy(std::begin(ref_u), std::end(ref_u), v.begin());
    }

    // manually calculate Bphi from the equivalent reference u
    B->Apply(*u, *y);

    // test y versus Bphi at l calculated from the fission sites
    {
        auto ref = Traits::get_data(y);
        auto v   = Traits::get_data(Bpsi_l);

        EXPECT_VEC_SOFTEQ(ref, v, 1.0e-5);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(FM_AccelerationTest, solve)
{
    using std::string;

    // setup the spn problem
    this->builder.setup("mesh4x4.xml");

    // build the spn problem
    this->acceleration->build_problem(this->builder);
    EXPECT_EQ(16, this->acceleration->mesh().num_cells());
    this->print_matrices();

    // initialize the spn problem
    this->acceleration->initialize(this->mc_db);

    // build fission source
    this->build_fs();

    // starting cycle
    this->acceleration->start_cycle(1.0, this->fs);

    // build l+1/2 fission source
    this->build_fs(true);

    // ending cycle
    this->acceleration->end_cycle(this->fs);
}

//---------------------------------------------------------------------------//
/*

  Row 0-1 on node 0
  Row 2-3 on node 1
  Row 4   on node 2
  Row 5   on node 3

       0   1   2   3   4   5

  0   10   2   0  31   2   3

  1    0   3  10   0   9   0

  2    9   5   1  11   3  13

  3    8   3   6  12   4  40

  4    1   1  12   0   2  14

  5    0   7   9  14   7   5


 */
TEST(Epetra, transpose)
{
#ifdef COMM_MPI
    if (profugus::nodes() != 4)
        return;
    int node = profugus::node();

    Epetra_MpiComm comm(profugus::communicator);

    std::vector<int> l2g;

    if (node == 0)
    {
        l2g = {0, 1};
    }
    if (node == 1)
    {
        l2g = {2, 3};
    }
    if (node == 2)
    {
        l2g = {4};
    }
    if (node == 3)
    {
        l2g = {5};
    }

    Epetra_Map map(-1, l2g.size(), &l2g[0], 0, comm);
    Epetra_CrsMatrix m(Copy, map, 0);
    Epetra_Vector x(map);
    {
        if (node == 0)
        {
            int i0[]    = {0, 1, 3, 4, 5};
            double v0[] = {10, 2, 31, 2, 3};

            m.InsertGlobalValues(0, 5, v0, i0);

            int i1[]    = {1, 2, 4};
            double v1[] = {3, 10, 9};

            m.InsertGlobalValues(1, 3, v1, i1);

            x[0] = 1.0;
            x[1] = 2.0;
        }

        if (node == 1)
        {
            int i0[]    = {0, 1, 2, 3, 4, 5};
            double v0[] = {9, 5, 1, 11, 3, 13};

            m.InsertGlobalValues(2, 6, v0, i0);

            int i1[]    = {0, 1, 2, 3, 4, 5};
            double v1[] = {8, 3, 6, 12, 4, 40};

            m.InsertGlobalValues(3, 6, v1, i1);

            x[0] = 3.0;
            x[1] = 4.0;
        }

        if (node == 2)
        {
            int i0[]    = {0, 1, 2, 4, 5};
            double v0[] = {1, 1, 12, 2, 14};

            m.InsertGlobalValues(4, 5, v0, i0);

            x[0] = 5.0;
        }

        if (node == 3)
        {
            int i0[]    = {1, 2, 3, 4, 5};
            double v0[] = {7, 9, 14, 7, 5};

            m.InsertGlobalValues(5, 5, v0, i0);

            x[0] = 6.0;
        }

        m.FillComplete();
    }

    // make transpose
    Epetra_CrsMatrix t(Copy, m.ColMap(), 0);
    {
        EXPECT_FALSE(t.IndicesAreLocal());

        // loop over local rows
        for (int row = 0; row < m.NumMyRows(); ++row)
        {
            int i[] = {m.GRID(row)};

            int n = m.NumMyCols();
            int l = 0;

            std::vector<int>    idx(n, 0);
            std::vector<double> val(n, 0.0);

            m.ExtractGlobalRowCopy(i[0], n, l, &val[0], &idx[0]);

            for (int col = 0; col < l; ++col)
            {
                int j = idx[col];

                double v[] = {val[col]};

                t.InsertGlobalValues(j, 1, v, i);
            }
        }

        t.FillComplete(m.OperatorRangeMap(), m.OperatorDomainMap());
    }

    // do regular mat-vec
    {
        Epetra_Vector y(map);
        m.Apply(x, y);

        if (node == 0)
        {
            EXPECT_EQ(166.0, y[0]);
            EXPECT_EQ(81.0, y[1]);
        }

        if (node == 1)
        {
            EXPECT_EQ(159.0, y[0]);
            EXPECT_EQ(340.0, y[1]);
        }

        if (node == 2)
        {
            EXPECT_EQ(133.0, y[0]);
        }

        if (node == 3)
        {
            EXPECT_EQ(162.0, y[0]);
        }
    }

    // do mat-vec with transpose option on m
    {
        Epetra_Vector y(map);
        m.SetUseTranspose(true);
        m.Apply(x, y);

        if (node == 0)
        {
            EXPECT_EQ(74.0, y[0]);
            EXPECT_EQ(82.0, y[1]);
        }

        if (node == 1)
        {
            EXPECT_EQ(161.0, y[0]);
            EXPECT_EQ(196.0, y[1]);
        }

        if (node == 2)
        {
            EXPECT_EQ(97.0, y[0]);
        }

        if (node == 3)
        {
            EXPECT_EQ(302.0, y[0]);
        }
    }

    // do mat-vec with transposed matrix
    {
        Epetra_Vector y(map);
        t.Apply(x, y);

        if (node == 0)
        {
            EXPECT_EQ(74.0, y[0]);
            EXPECT_EQ(82.0, y[1]);
        }

        if (node == 1)
        {
            EXPECT_EQ(161.0, y[0]);
            EXPECT_EQ(196.0, y[1]);
        }

        if (node == 2)
        {
            EXPECT_EQ(97.0, y[0]);
        }

        if (node == 3)
        {
            EXPECT_EQ(302.0, y[0]);
        }
    }

    // Transpose with the transposer
    Epetra_RowMatrixTransposer trans(&m);
    Epetra_CrsMatrix *p;

    trans.CreateTranspose(true, p);
    {
        Epetra_Vector y(map);
        p->Apply(x, y);

        if (node == 0)
        {
            EXPECT_EQ(74.0, y[0]);
            EXPECT_EQ(82.0, y[1]);
        }

        if (node == 1)
        {
            EXPECT_EQ(161.0, y[0]);
            EXPECT_EQ(196.0, y[1]);
        }

        if (node == 2)
        {
            EXPECT_EQ(97.0, y[0]);
        }

        if (node == 3)
        {
            EXPECT_EQ(302.0, y[0]);
        }
    }

    delete p;

#endif
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Matrix_Acceleration.cc
//---------------------------------------------------------------------------//
