//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstShiftedInverseOperator.cc
 * \author Steven Hamilton
 * \brief  ShiftedInverseOperator test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <SPn/config.h>

#include "harness/DBC.hh"
#include "../ShiftedInverseOperator.hh"

#include "LinAlgTraits.hh"

using namespace std;

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class ShiftedInverseTest: public ::testing::Test
{
  protected:

    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_ParameterList;
    typedef Epetra_MultiVector                      MV;
    typedef Epetra_Operator                         OP;
    typedef Epetra_CrsMatrix                        Matrix;
    typedef profugus::ShiftedInverseOperator<MV,OP> ShiftInvOp;

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
        int my_size = global_size / nodes;

        // Build CrsMatrix
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",global_size);
        d_B = linalg_traits::build_matrix<Matrix>("diagonal",global_size);

        // Create options database
        d_db = Teuchos::rcp(new Teuchos::ParameterList("test"));
        d_db->set("solver_type", std::string("stratimikos"));
        d_db->set("tolerance",1e-10);
        d_db->set("max_itr",20);

        // Build solver
        d_operator = Teuchos::rcp(new ShiftInvOp(d_db));
        CHECK(!d_operator.is_null());
    }

  protected:
    int node;
    int nodes;

    RCP_ParameterList              d_db;
    Teuchos::RCP<Epetra_CrsMatrix> d_A;
    Teuchos::RCP<Epetra_CrsMatrix> d_B;
    Teuchos::RCP<ShiftInvOp>       d_operator;
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEST_F(ShiftedInverseTest, basic)
{
    Teuchos::RCP<MV> x1 = linalg_traits::build_vector<MV>(d_A->NumGlobalRows());
    Teuchos::RCP<MV> x2 = linalg_traits::build_vector<MV>(d_A->NumGlobalRows());
    Teuchos::RCP<MV>  y = linalg_traits::build_vector<MV>(d_A->NumGlobalRows());

    // Solver tolerance is 1e-10, this should give is 1e-8 in vector entries
    double tol = 1e-8;

    // Test all cases against two vectors
    x1->PutScalar(1.0);
    for( int i=0; i<x2->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        (*x2)[0][i] = static_cast<double>(global_row+1);
    }

    // First set only one operator
    d_operator->set_operator(d_A);
    d_operator->set_shift(0.0);

    // First vector
    d_operator->Apply(*x1,*y);

    double ref[] = {
        1.0000e+01, 1.9000e+01, 2.7000e+01, 3.4000e+01, 4.0000e+01,
        4.5000e+01, 4.9000e+01, 5.2000e+01, 5.4000e+01, 5.5000e+01,
        5.5000e+01, 5.4000e+01, 5.2000e+01, 4.9000e+01, 4.5000e+01,
        4.0000e+01, 3.4000e+01, 2.7000e+01, 1.9000e+01, 1.0000e+01 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref[global_row], (*y)[0][i], tol );
    }

    // Second vector
    d_operator->Apply(*x2,*y);

    double ref2[] = {
        7.33333333e+01, 1.45666667e+02, 2.16000000e+02, 2.83333333e+02,
        3.46666667e+02, 4.05000000e+02, 4.57333333e+02, 5.02666667e+02,
        5.40000000e+02, 5.68333333e+02, 5.86666667e+02, 5.94000000e+02,
        5.89333333e+02, 5.71666667e+02, 5.40000000e+02, 4.93333333e+02,
        4.30666667e+02, 3.51000000e+02, 2.53333333e+02, 1.36666667e+02 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref2[global_row], (*y)[0][i], tol );
    }

    // Set nonzero shift
    d_operator->set_shift(0.5);

    // First vector
    d_operator->Apply(*x1,*y);

    double ref3[] = {
        4.36933798e+00,  5.55400697e+00,  2.96167247e+00, -2.11149826e+00,
       -7.12891986e+00, -9.58188153e+00, -8.24390244e+00, -3.78397213e+00,
        1.56794425e+00,  5.13588850e+00,  5.13588850e+00,  1.56794425e+00,
       -3.78397213e+00, -8.24390244e+00, -9.58188153e+00, -7.12891986e+00,
       -2.11149826e+00,  2.96167247e+00,  5.55400697e+00,  4.36933798e+00 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref3[global_row], (*y)[0][i], tol );
    }

    // Second vector
    d_operator->Apply(*x2,*y);

    double ref4[] = {
        5.29016624e+01,  7.83524936e+01,  6.26270780e+01,  1.25881234e+01,
       -4.77448929e+01, -8.92054627e+01, -9.20633012e+01, -5.58894891e+01,
        2.29067586e-01,  4.72330904e+01,  6.06205681e+01,  3.26977617e+01,
       -2.35739256e+01, -8.10586500e+01, -1.12014049e+02, -1.01962424e+02,
       -5.69295868e+01, -4.31956019e-01,  3.82816528e+01,  3.88544352e+01 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref4[global_row], (*y)[0][i], tol );
    }

    // Set second operator and shift=0, results should be same as original
    d_operator->set_rhs_operator(d_B);
    d_operator->set_shift(0.0);

    // First vector
    d_operator->Apply(*x1,*y);

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref[global_row], (*y)[0][i], tol );
    }

    // Second vector
    d_operator->Apply(*x2,*y);

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref2[global_row], (*y)[0][i], tol );
    }

    // Nonzero shift
    d_operator->set_shift(0.5);

    // First vector
    d_operator->Apply(*x1,*y);

    double ref5[] = {
        1.84013100e+00,  1.76019650e+00, -1.07993450e+00, -3.30016375e+00,
        7.99345003e-02,  2.26019650e+00, -3.34013100e+00,  1.75000000e+00,
       -1.15986900e+00,  1.49672501e-01, -2.89148504e-01, -1.37652738e-01,
       -1.60240543e-01, -1.41264818e-01, -1.33435369e-01, -1.24840654e-01,
       -1.17520705e-01, -1.11274761e-01, -1.03555965e-01, -1.12055504e-01 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref5[global_row], (*y)[0][i], tol );
    }

    // Second vector
    d_operator->Apply(*x2,*y);

    double ref6[] = {
        1.06031307e+00,  5.90469609e-01, -2.46984346e+00, -4.82539134e+00,
       -1.53015654e+00,  5.90469609e-01, -5.06031307e+00,  0.00000000e+00,
       -2.93968693e+00, -1.65078268e+00, -2.10796503e+00, -1.97133972e+00,
       -2.00667608e+00, -1.99861790e+00, -2.00023440e+00, -2.00009290e+00,
       -1.99920823e+00, -2.00505362e+00, -1.96541646e+00, -2.25432294e+00 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref6[global_row], (*y)[0][i], tol );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstShiftedInverseOperator.cc
//---------------------------------------------------------------------------//
