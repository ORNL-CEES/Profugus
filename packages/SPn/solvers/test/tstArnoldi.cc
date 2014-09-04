//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstArnoldi.cc
 * \author Steven Hamilton
 * \date   Mon Aug 22 12:22:41 2011
 * \brief  Test Arnoldi class.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>

#include <SPn/config.h>

#include "../Arnoldi.hh"
#include "../InverseOperator.hh"
#include "LinAlgTraits.hh"

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

using namespace std;

//---------------------------------------------------------------------------//
// Test Fixture
//---------------------------------------------------------------------------//

class Arnoldi_Test : public testing::Test
{
  protected:

#ifdef COMM_MPI
    typedef Epetra_MpiComm      Comm;
#else
    typedef Epetra_SerialComm   Comm;
#endif

    typedef Epetra_MultiVector               MV;
    typedef Epetra_Operator                  OP;
    typedef profugus::Arnoldi<MV,OP>         Arnoldi;
    typedef Epetra_Map                       Map;
    typedef Epetra_CrsGraph                  Graph;
    typedef Epetra_CrsMatrix                 Matrix;
    typedef Arnoldi::RCP_MV                  RCP_MV;
    typedef Arnoldi::RCP_OP                  RCP_OP;
    typedef Teuchos::RCP<Matrix>             RCP_Matrix;
    typedef profugus::InverseOperator<MV,OP> InverseOperator;
    typedef Arnoldi::RCP_ParameterList       RCP_ParameterList;
    typedef Arnoldi::ParameterList           ParameterList;
    typedef std::vector<int>                 Vec_Int;
    typedef std::vector<double>              Vec_Dbl;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        // Create Epetra MPI Comm(he's not the first person to observe this)
#ifdef COMM_MPI
        comm = Teuchos::rcp(new Comm(MPI_COMM_WORLD));
#else
        comm = Teuchos::rcp(new Comm);
#endif
        db  = Teuchos::rcp(new ParameterList("test"));
        db2 = Teuchos::rcp(new ParameterList("test"));
    }

  protected:
    int node, nodes;
    Teuchos::RCP<Comm> comm;
    RCP_ParameterList  db, db2;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Arnoldi_Test, Eigensolver)
{
    using Teuchos::RCP;

    if (nodes > 2)
        return;

    // Matrix size
    int N = 8;

    RCP_Matrix A = linalg_traits::build_matrix<Matrix>("laplacian",N);
    RCP_Matrix B = linalg_traits::build_matrix<Matrix>("diagonal",N);
    RCP_MV e_vec = linalg_traits::build_vector<MV>(N);
    e_vec->Random();

    // Test with matrix A
    double e_val;

    // database settings
    {
        // make anasazi settings
        ParameterList adb("anasazi");
        adb.set("Convergence Tolerance", 1.0e-10);
        adb.set("Maximum Restarts", 50);
        db->set("Anasazi", adb);
    }

    // make solver
    Arnoldi solver( db );
    solver.set_operator(A);
    // This is actually not a great test problem because the dominant
    // eigenvector is exactly orthogonal to a constant vector.  We have
    // to initialize with a random vector or else we'll converge to the
    // 2nd eigenmode.
    solver.solve( e_val, e_vec );

    cout.precision(15);
    cout << "Eig(A) = " << e_val << endl;

    int offset;
    if( node==0 )
        offset=0;
    else
        offset=N/2;

    for( int i=0; i<e_vec->MyLength(); ++i )
        cout << "Evec[" << i+offset << "]: " << (*e_vec)[0][i] << endl;

    vector<double> eref(8);
    eref[0] =   0.161229841765317;
    eref[1] =  -0.303012985114695;
    eref[2] =   0.408248290463863;
    eref[3] =  -0.464242826880012;
    eref[4] =   0.464242826880013;
    eref[5] =  -0.408248290463863;
    eref[6] =   0.303012985114696;
    eref[7] =  -0.161229841765317;

    // Test against Matlab computed values
    double tol = 1.0e-6;
    EXPECT_SOFTEQ(e_val, 3.879385241571816, tol);

    // Get sign of computed eigenvector
    double sign = 1.0;
    if( (*e_vec)[0][0] < 0.0 )
        sign = -1.0;

    // Check eigenvector values
    // Can only check absolute values because e-vec can vary +/-
    //  (running with a different number of processors may change the sign)
    for( int i=0; i<e_vec->MyLength(); ++i )
        EXPECT_SOFTEQ(sign*(*e_vec)[0][i], eref[i+offset], tol);

    // Create operator for A^{-1}B
    {
        db2->set("max_itr", 50);
        db2->set("tolerance", 1e-12);
        db2->set("solver_type", std::string("stratimikos"));
    }

    RCP<InverseOperator> AinvB(Teuchos::rcp(new InverseOperator(db2)));
    AinvB->set_operator(A);
    AinvB->set_rhs_operator(B);

    // Now an eigensolver for A^{-1}B
    Arnoldi gensolver( db2 );
    gensolver.set_operator( AinvB );

    // Solve
    e_vec = linalg_traits::build_vector<MV>(N);
    gensolver.solve( e_val, e_vec );

    cout << "Eig(AinvB) = " << e_val << endl;

    for( int i=0; i<e_vec->MyLength(); ++i )
        cout << "Evec[" << i+offset << "]: " << (*e_vec)[0][i] << endl;

    eref[0] = 0.125730096111867;
    eref[1] = 0.248228875407670;
    eref[2] = 0.357968479865590;
    eref[3] = 0.440108258878819;
    eref[4] = 0.477004161945638;
    eref[5] = 0.452604018262593;
    eref[6] = 0.358411182752022;
    eref[7] = 0.199739110973060;

    // Test against Matlab computed values
    EXPECT_SOFTEQ(e_val, 38.909863460868493, tol);

    // Get sign of computed eigenvector
    sign = 1.0;
    if( (*e_vec)[0][0] < 0.0 )
        sign = -1.0;

    // Check eigenvector values
    for( int i=0; i<e_vec->MyLength(); ++i )
    {
        EXPECT_SOFTEQ(sign*(*e_vec)[0][i], eref[i+offset], tol);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstArnoldi.cc
//---------------------------------------------------------------------------//
