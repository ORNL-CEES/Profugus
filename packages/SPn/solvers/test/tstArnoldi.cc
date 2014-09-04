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
#include <random>

#include <SPn/config.h>

#include "../Arnoldi.hh"
#include "../InverseOperator.hh"
#include "LinAlgTraits.hh"

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

using namespace std;

//---------------------------------------------------------------------------//
// Test Fixture
//---------------------------------------------------------------------------//

class Arnoldi_Test : public testing::Test
{
  protected:

    typedef Epetra_MultiVector               MV;
    typedef Epetra_Operator                  OP;
    typedef Epetra_CrsMatrix                 Matrix;
    typedef profugus::Arnoldi<MV,OP>         Arnoldi;
    typedef profugus::InverseOperator<MV,OP> InverseOperator;
    typedef Arnoldi::RCP_ParameterList       RCP_ParameterList;
    typedef Arnoldi::ParameterList           ParameterList;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        db  = Teuchos::rcp(new ParameterList("test"));
        db2 = Teuchos::rcp(new ParameterList("test"));
    }

  protected:
    int node, nodes;
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

    Teuchos::RCP<Matrix> A = linalg_traits::build_matrix<Matrix>("laplacian",N);
    Teuchos::RCP<Matrix> B = linalg_traits::build_matrix<Matrix>("diagonal",N);
    Teuchos::RCP<MV> e_vec = linalg_traits::build_vector<MV>(N);
    random_device rd;
    std::vector<double> init(N);
    for( int i=0; i<N; ++i )
        init[i] = rd();
    linalg_traits::fill_vector<MV>(e_vec,init);

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

    vector<double> eref =
      { 0.161229841765317, -0.303012985114695, 0.408248290463863,
       -0.464242826880012, 0.464242826880013, -0.408248290463863,
        0.303012985114696, -0.161229841765317};

    // Test against Matlab computed values
    double tol = 1.0e-6;
    EXPECT_SOFTEQ(e_val, 3.879385241571816, tol);

    linalg_traits::set_sign(e_vec);
    linalg_traits::test_vector<MV>(e_vec,eref);

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
    Arnoldi gensolver( db );
    gensolver.set_operator( AinvB );

    // Solve
    e_vec = linalg_traits::build_vector<MV>(N);
    gensolver.solve( e_val, e_vec );

    cout << "Eig(AinvB) = " << e_val << endl;

    std::vector<double> eref2 =
        {0.125730096111867, 0.248228875407670, 0.357968479865590,
         0.440108258878819, 0.477004161945638, 0.452604018262593,
         0.358411182752022, 0.199739110973060};

    // Test against Matlab computed values
    EXPECT_SOFTEQ(e_val, 38.909863460868493, tol);

    linalg_traits::set_sign(e_vec);
    linalg_traits::test_vector<MV>(e_vec,eref2);
}

//---------------------------------------------------------------------------//
//                        end of tstArnoldi.cc
//---------------------------------------------------------------------------//
