//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstPowerIteration.cc
 * \author Steven Hamilton
 * \brief  PowerIteration unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <SPn/config.h>

#include "../PowerIteration.hh"

#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class PowerIterationTest : public ::testing::Test
{
  protected:

    typedef Epetra_MultiVector                MV;
    typedef Epetra_Operator                   OP;
    typedef Epetra_CrsMatrix                  Matrix;
    typedef profugus::PowerIteration<MV,OP>   PowerIteration;
    typedef PowerIteration::RCP_ParameterList RCP_ParameterList;
    typedef PowerIteration::ParameterList     ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Parallelism
        node  = profugus::node();
        nodes = profugus::nodes();

        // Build an Epetra map
        d_N = 20;
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",d_N);

        // Build eigenvector
        d_x = linalg_traits::build_vector<MV>(d_N);
        std::vector<double> one(d_N,1.0);
        linalg_traits::fill_vector<MV>(d_x,one);

        // Create options database
        d_db = Teuchos::rcp(new ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);

        // Build solver
        d_solver = Teuchos::rcp(new PowerIteration(d_db));
        CHECK(!d_solver.is_null());
        d_solver->set_operator(d_A);
    }

    void solve()
    {
        d_solver->solve(d_lambda,d_x);
        d_iters = d_solver->num_iters();
        d_converged = d_solver->converged();
    }

  protected:
    int node;
    int nodes;
    int d_N;

    RCP_ParameterList            d_db;
    Teuchos::RCP<Matrix>         d_A;
    Teuchos::RCP<MV>             d_x;
    Teuchos::RCP<PowerIteration> d_solver;
    double                       d_lambda;

    int d_iters;
    bool d_converged;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST_F(PowerIterationTest, basic)
{
    // Run two iterations and stop
    solve();

    // Make sure solver reports that two iterations were performed
    //  and that it is not converged
    EXPECT_EQ( 2, d_iters );
    EXPECT_TRUE( !d_converged );

    // Reset initial vector and re-solve
    d_solver->set_max_iters(1000);
    std::vector<double> one(d_N,1.0);
    linalg_traits::fill_vector<MV>(d_x,one);
    solve();

    EXPECT_EQ( 261, d_iters );
    EXPECT_TRUE( d_converged );
    EXPECT_SOFTEQ( d_lambda, 3.911145611572282, 1.0e-6 );

    // Compare against reference solution from Matlab
    std::vector<double> ref =
            {9.096342209328087e-02,  -1.738443448352310e-01,
             2.412784344502512e-01,  -2.872738761376646e-01,
             3.077437731144456e-01,  -3.008692853090870e-01,
             2.672612412471236e-01,  -2.099058632158057e-01,
             1.338994282792137e-01,  -4.599544168741217e-02,
             -4.599544168741519e-02,  1.338994282792167e-01,
             -2.099058632158085e-01,  2.672612412471261e-01,
             -3.008692853090894e-01,  3.077437731144476e-01,
             -2.872738761376664e-01,  2.412784344502525e-01,
             -1.738443448352319e-01,  9.096342209328130e-02};

    linalg_traits::test_vector<MV>(d_x,ref);

    // Solve again, should return without iterating
    solve();

    EXPECT_EQ( 1, d_iters );
    EXPECT_TRUE( d_converged );

    // Make sure solution didn't change
    linalg_traits::test_vector<MV>(d_x,ref);
}

//---------------------------------------------------------------------------//
//                        end of tstPowerIteration.cc
//---------------------------------------------------------------------------//
