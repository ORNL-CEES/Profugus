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

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class RichardsonTest : public testing::Test
{
  protected:

    typedef Epetra_MultiVector            MV;
    typedef Epetra_Operator               OP;
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

        // Build Epetra communicator
#ifdef COMM_MPI
        Epetra_MpiComm comm(profugus::communicator);
#else
        Epetra_SerialComm comm;
#endif

        // Build an Epetra map
        int global_size = 8;
        d_map = Teuchos::rcp( new Epetra_Map( global_size, 0, comm ) );

        int my_size = global_size / nodes;

        std::cout << "Global size: " << global_size << std::endl;
        std::cout << "Local size: " << my_size << std::endl;

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
                vals[0] =  0.50;
                vals[1] = -0.25;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else if( global_row == global_size-1 )
            {
                std::vector<int> ids(2);
                ids[0] = global_size-2;
                ids[1] = global_size-1;
                std::vector<double> vals(2);
                vals[0] = -0.25;
                vals[1] =  0.50;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else
            {
                std::vector<int> ids(3);
                ids[0] = global_row-1;
                ids[1] = global_row;
                ids[2] = global_row+1;
                std::vector<double> vals(3);
                vals[0] = -0.25;
                vals[1] =  0.50;
                vals[2] = -0.25;
                d_A->InsertGlobalValues(global_row,3,&vals[0],&ids[0]);
            }
        }
        d_A->FillComplete();

        std::cout << "Matrix built" << std::endl;

        // Build lhs and rhs vectors
        d_x = Teuchos::rcp( new Epetra_MultiVector(*d_map,1) );
        d_b = Teuchos::rcp( new Epetra_MultiVector(*d_map,1) );
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_map->GID(my_row);
            d_b->ReplaceGlobalValue(global_row,0,
                    static_cast<double>(8-global_row));
        }

        // Create options database
        d_db = Teuchos::rcp(new ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);
        d_db->set("Damping Factor",2.0);

        // Build solver
        d_solver = Teuchos::rcp(new Richardson(d_db));
        Check(!d_solver.is_null());
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
    Teuchos::RCP<Epetra_Map>         d_map;
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
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref[global_row], (*d_x)[0][my_row], 1.0e-7 );
    }

    // Solve again, should return without iterating
    solve();

    EXPECT_EQ( 0, d_iters );
    EXPECT_TRUE( d_converged );

    // Make sure solution didn't change
    for( int my_row = 0; my_row < 10/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref[global_row], (*d_x)[0][my_row], 1.0e-7 );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstRichardson.cc
//---------------------------------------------------------------------------//
