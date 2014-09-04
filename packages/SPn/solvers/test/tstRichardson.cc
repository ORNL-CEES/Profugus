//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstRichardson.cc
 * \author Steven Hamilton
 * \brief  Richardson unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include <SPn/config.h>
#include "../Richardson.hh"

#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class RichardsonTest : public testing::Test
{
  protected:

    typedef Epetra_MultiVector            MV;
    typedef Epetra_Operator               OP;
    typedef Epetra_CrsMatrix              Matrix;
    typedef profugus::Richardson<MV,OP>   Richardson;
    typedef Richardson::RCP_ParameterList RCP_ParameterList;
    typedef Richardson::ParameterList     ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Parallelism
        node  = profugus::node();
        nodes = profugus::nodes();

        // Build an Epetra map
        int global_size = 8;
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",global_size);
        int my_size = global_size / nodes;

        // Build lhs and rhs vectors
        d_x = linalg_traits::build_vector<MV>(global_size);
        d_b = linalg_traits::build_vector<MV>(global_size);
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_A->GRID(my_row);
            d_b->ReplaceGlobalValue(global_row,0,
                    static_cast<double>(4*(8-global_row)));
        }

        // Create options database
        d_db = Teuchos::rcp(new ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);
        d_db->set("Damping Factor",0.5);

        // Build solver
        d_solver = Teuchos::rcp(new Richardson(d_db));
        CHECK(!d_solver.is_null());
        d_solver->set_operator(d_A);
    }

    void solve()
    {
        d_solver->solve(d_x,d_b);
        d_iters = d_solver->num_iters();
        d_converged = d_solver->converged();
    }

  protected:
    int node;
    int nodes;

    RCP_ParameterList                d_db;
    Teuchos::RCP<Epetra_CrsMatrix>   d_A;
    Teuchos::RCP<Epetra_MultiVector> d_x;
    Teuchos::RCP<Epetra_MultiVector> d_b;
    Teuchos::RCP<Richardson>         d_solver;

    int d_iters;
    bool d_converged;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RichardsonTest, basic)
{
    // Run two iterations and stop
    solve();

    // Make sure solver reports that two iterations were performed
    //  and that it is not converged
    EXPECT_EQ( 2, d_iters );
    EXPECT_TRUE( !d_converged );

    // Reset initial vector and re-solve
    d_solver->set_max_iters(1000);
    d_x->PutScalar(0.0);
    solve();

    EXPECT_EQ( 294, d_iters ); // Iteration count from Matlab implementation
    EXPECT_TRUE( d_converged );

    // Compare against reference solution from Matlab
    double ref[] = {
        90.6666666666667,
       149.3333333333333,
       180.0000000000000,
       186.6666666666666,
       173.3333333333333,
       144.0000000000000,
       102.6666666666666,
        53.3333333333333};

    for( int my_row = 0; my_row < 10/nodes; ++my_row )
    {
        int global_row = d_A->GRID(my_row);
        EXPECT_SOFTEQ( ref[global_row], (*d_x)[0][my_row], 1.0e-7 );
    }

    // Solve again, should return without iterating
    solve();

    EXPECT_EQ( 0, d_iters );
    EXPECT_TRUE( d_converged );

    // Make sure solution didn't change
    for( int my_row = 0; my_row < 10/nodes; ++my_row )
    {
        int global_row = d_A->GRID(my_row);
        EXPECT_SOFTEQ( ref[global_row], (*d_x)[0][my_row], 1.0e-7 );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstRichardson.cc
//---------------------------------------------------------------------------//
