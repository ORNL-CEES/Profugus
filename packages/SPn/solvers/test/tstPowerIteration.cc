//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstPowerIteration.cc
 * \author Steven Hamilton
 * \brief  PowerIteration unit test.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <SPn/config.h>

#include "../PowerIteration.hh"

#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class PowerIterationTest : public ::testing::Test
{
  protected:

    typedef Epetra_MultiVector                MV;
    typedef Epetra_Operator                   OP;
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

        // Build Epetra communicator
#ifdef COMM_MPI
        Epetra_MpiComm comm(profugus::communicator);
#else
        Epetra_SerialComm comm;
#endif

        // Build an Epetra map
        int global_size = 20;
        d_map = Teuchos::rcp( new Epetra_Map( global_size, 0, comm ) );

        int my_size = global_size / nodes;

        // Build CrsMatrix
        d_A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,*d_map,3) );
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_map->GID(my_row);
            if( global_row == 0 )
            {
                std::vector<int> ids(2);
                ids[0] = 0;
                ids[1] = 1;
                std::vector<double> vals(2);
                vals[0] =  2.0;
                vals[1] = -1.0;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else if( global_row == global_size-1 )
            {
                std::vector<int> ids(2);
                ids[0] = 18;
                ids[1] = 19;
                std::vector<double> vals(2);
                vals[0] = -1.0;
                vals[1] =  2.0;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else
            {
                std::vector<int> ids(3);
                ids[0] = global_row-1;
                ids[1] = global_row;
                ids[2] = global_row+1;
                std::vector<double> vals(3);
                vals[0] = -1.0;
                vals[1] =  2.0;
                vals[2] = -1.0;
                d_A->InsertGlobalValues(global_row,3,&vals[0],&ids[0]);
            }
        }
        d_A->FillComplete();

        // Build eigenvector
        d_x = Teuchos::rcp( new Epetra_MultiVector(*d_map,1) );
        d_x->PutScalar(1.0);

        // Create options database
        d_db = Teuchos::rcp(new ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);

        // Build solver
        d_solver = Teuchos::rcp(new PowerIteration(d_db));
        Check (!d_solver.is_null());
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

    RCP_ParameterList                d_db;
    Teuchos::RCP<Epetra_Map>         d_map;
    Teuchos::RCP<Epetra_CrsMatrix>   d_A;
    Teuchos::RCP<Epetra_MultiVector> d_x;
    Teuchos::RCP<Epetra_MultiVector> d_b;
    Teuchos::RCP<PowerIteration>     d_solver;
    double                           d_lambda;

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
    d_x->PutScalar(1.0);
    solve();

    EXPECT_EQ( 261, d_iters );
    EXPECT_TRUE( d_converged );
    EXPECT_SOFTEQ( d_lambda, 3.911145611572282, 1.0e-6 );

    // Compare against reference solution from Matlab
    double ref[] = {
             9.096342209328087e-02,  -1.738443448352310e-01,
             2.412784344502512e-01,  -2.872738761376646e-01,
             3.077437731144456e-01,  -3.008692853090870e-01,
             2.672612412471236e-01,  -2.099058632158057e-01,
             1.338994282792137e-01,  -4.599544168741217e-02,
             -4.599544168741519e-02,  1.338994282792167e-01,
             -2.099058632158085e-01,  2.672612412471261e-01,
             -3.008692853090894e-01,  3.077437731144476e-01,
             -2.872738761376664e-01,  2.412784344502525e-01,
             -1.738443448352319e-01,  9.096342209328130e-02};

    for( int my_row = 0; my_row < 20/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref[global_row], (*d_x)[0][my_row], 1.0e-6 );
    }

    // Solve again, should return without iterating
    solve();

    EXPECT_EQ( 1, d_iters );
    EXPECT_TRUE( d_converged );

    // Make sure solution didn't change
    for( int my_row = 0; my_row < 20/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref[global_row], (*d_x)[0][my_row], 1.0e-6 );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstPowerIteration.cc
//---------------------------------------------------------------------------//
