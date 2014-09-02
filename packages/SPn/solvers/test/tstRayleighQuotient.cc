//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstRayleighQuotient.cc
 * \author Steven Hamilton
 * \brief  RayleighQuotient unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>
#include <string>

#include <SPn/config.h>

#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include "../RayleighQuotient.hh"
#include "../ShiftedInverseOperator.hh"

// Reference solution from matlab

double ref_eigenvalue = 0.4890748754542557;
double ref_eigenvector[] = {
    4.599544191e-02, 9.096342166e-02, 1.338994289e-01, 1.738443441e-01,
    2.099058640e-01, 2.412784337e-01, 2.672612419e-01, 2.872738756e-01,
    3.008692856e-01, 3.077437730e-01, 3.077437730e-01, 3.008692856e-01,
    2.872738756e-01, 2.672612419e-01, 2.412784337e-01, 2.099058640e-01,
    1.738443441e-01, 1.338994289e-01, 9.096342166e-02, 4.599544191e-02 };

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class RQITest : public ::testing::Test
{
  protected:

    typedef Epetra_MultiVector                  MV;
    typedef Epetra_Operator                     OP;
    typedef profugus::RayleighQuotient<MV,OP>   RayleighQuotient;
    typedef RayleighQuotient::RCP_ParameterList RCP_ParameterList;
    typedef RayleighQuotient::ParameterList     ParameterList;

  protected:

    // Initialization that are performed for each test
    void build_solver()
    {
        using Teuchos::rcp;

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
        d_B = rcp( new Epetra_CrsMatrix(Copy,*d_map,1) );
        d_A = rcp( new Epetra_CrsMatrix(Copy,*d_map,3) );
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_map->GID(my_row);
            if( global_row == 0 )
            {
                std::vector<int> ids(2);
                ids[0] = 0;
                ids[1] = 1;
                std::vector<double> vals(2);
                vals[0] =  3.0;
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
                vals[1] =  3.0;
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
                vals[1] =  3.0;
                vals[2] = -1.0;
                d_A->InsertGlobalValues(global_row,3,&vals[0],&ids[0]);
            }
            std::vector<int>    inds(1);
            std::vector<double> vals(1);
            inds[0] = global_row;
            vals[0] = 0.5;
            d_B->InsertGlobalValues(global_row,1,&vals[0],&inds[0]);
        }
        d_A->FillComplete();
        d_B->FillComplete();

        // Build eigenvector
        d_x = Teuchos::rcp( new Epetra_MultiVector(*d_map,1) );

        // Create options database
        d_db = rcp(new ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);
        if( d_use_fixed_shift )
        {
            d_db->set("use_fixed_shift",d_use_fixed_shift);
            d_db->set("eig_shift",d_eig_shift);
        }
        RCP_ParameterList op_db = Teuchos::sublist(d_db, "operator_db");
        op_db->set("solver_type", std::string("stratimikos"));
        op_db->set("tolerance",1e-8);
        op_db->set("max_itr",20);

        // Create ShiftedInverseOperator
        Teuchos::RCP<profugus::ShiftedInverseOperator<MV,OP> > shift_op =
            rcp(new profugus::ShiftedInverseOperator<MV,OP>(op_db));
        CHECK(!shift_op.is_null());
        shift_op->set_operator(d_A);
        shift_op->set_rhs_operator(d_B);

        // Build solver
        d_solver = rcp(new RayleighQuotient(d_db));
        CHECK(!d_solver.is_null());
        d_solver->set_operator(d_A);
        d_solver->set_rhs_operator(d_B);
        d_solver->set_shifted_operator(shift_op);
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
    Teuchos::RCP<Epetra_CrsMatrix>   d_B;
    Teuchos::RCP<Epetra_MultiVector> d_x;
    Teuchos::RCP<RayleighQuotient>   d_solver;
    double                           d_lambda;

    int d_iters;
    bool d_converged;
    bool d_use_fixed_shift;
    double d_eig_shift;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST_F(RQITest, basic)
{
    double eig_tol = 1e-10;
    double vec_tol = 1e-7;
    d_use_fixed_shift = false;
    build_solver();

    // Run two iterations and stop
    d_x->PutScalar(1.0);
    d_lambda = 1.0;
    solve();

    // Make sure solver reports that two iterations were performed
    //  and that it is not converged
    EXPECT_EQ( 2, d_iters );
    EXPECT_TRUE( !d_converged );

    // Reset initial vector and re-solve
    d_solver->set_max_iters(10);
    d_x->PutScalar(1.0);
    d_lambda = 1.0;
    solve();

    EXPECT_EQ( 4, d_iters );
    EXPECT_TRUE( d_converged );
    EXPECT_SOFTEQ( d_lambda, ref_eigenvalue, eig_tol );

    double sign = (*d_x)[0][0] / std::fabs((*d_x)[0][0]);
    for( int my_row = 0; my_row < 20/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref_eigenvector[global_row],
                       sign*(*d_x)[0][my_row], vec_tol );
    }

    // Solve again, should return in 1 iteration
    solve();

    EXPECT_EQ( 1, d_iters );
    EXPECT_TRUE( d_converged );
    EXPECT_SOFTEQ( d_lambda, ref_eigenvalue, eig_tol );

    // Make sure solution didn't change
    sign = (*d_x)[0][0] / std::fabs((*d_x)[0][0]);
    for( int my_row = 0; my_row < 20/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref_eigenvector[global_row],
                       sign*(*d_x)[0][my_row], vec_tol );
    }

    // Now reset and solve with fixed shift
    d_use_fixed_shift = true;
    d_eig_shift = 0.5;
    build_solver();

    // Run two iterations and stop
    d_x->PutScalar(1.0);
    d_lambda = 1.0;
    solve();

    // Make sure solver reports that two iterations were performed
    //  and that it is not converged
    EXPECT_EQ( 2, d_iters );
    EXPECT_TRUE( !d_converged );

    // Reset initial vector and re-solve
    d_solver->set_max_iters(1000);
    d_x->PutScalar(1.0);
    d_lambda = 1.0;
    solve();

    EXPECT_EQ( 8, d_iters ); // Heuristic
    EXPECT_TRUE( d_converged );
    EXPECT_SOFTEQ( d_lambda, ref_eigenvalue, eig_tol );

    sign = (*d_x)[0][0] / std::fabs((*d_x)[0][0]);
    for( int my_row = 0; my_row < 20/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref_eigenvector[global_row],
                       sign*(*d_x)[0][my_row], vec_tol );
    }

    // Solve again, should return in 1 iteration
    solve();

    EXPECT_EQ( 1, d_iters );
    EXPECT_TRUE( d_converged );
    EXPECT_SOFTEQ( d_lambda, ref_eigenvalue, eig_tol );

    // Make sure solution didn't change
    sign = (*d_x)[0][0] / std::fabs((*d_x)[0][0]);
    for( int my_row = 0; my_row < 20/nodes; ++my_row )
    {
        int global_row = d_map->GID(my_row);
        EXPECT_SOFTEQ( ref_eigenvector[global_row],
                       sign*(*d_x)[0][my_row], vec_tol );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstRayleighQuotient.cc
//---------------------------------------------------------------------------//
