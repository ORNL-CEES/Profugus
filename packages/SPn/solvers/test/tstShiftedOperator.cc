//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstShiftedOperator.cc
 * \author Steven Hamilton
 * \brief  ShiftedOperator unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../ShiftedOperator.hh"

#include <SPn/config.h>

#include "comm/global.hh"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class ShiftedOperatorTest : public testing::Test
{
  protected:
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    typedef Epetra_CrsMatrix   Matrix;
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

        // Build eigenvector
        d_x = linalg_traits::build_vector<MV>(global_size);
        d_y = linalg_traits::build_vector<MV>(global_size);

        // Build solver
        d_operator = Teuchos::rcp(new profugus::ShiftedOperator<MV,OP>());
        CHECK(!d_operator.is_null());
        d_operator->set_operator(d_A);
        d_operator->set_rhs_operator(d_B);
    }

  protected:
    int node;
    int nodes;

    Teuchos::RCP<Matrix>                  d_A;
    Teuchos::RCP<Matrix>                  d_B;
    Teuchos::RCP<MV>                      d_x;
    Teuchos::RCP<MV>                      d_y;
    Teuchos::RCP<profugus::ShiftedOperator<MV,OP> > d_operator;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ShiftedOperatorTest, basic)
{
    // Unshifted operator, same as multiplying by A
    d_operator->set_shift(0.0);
    d_x->PutScalar(1.0);
    d_operator->Apply(*d_x,*d_y);

    for( int i=0; i<d_y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        if( global_row == 0 || global_row == 19 )
        {
            EXPECT_SOFTEQ( 1.0, (*d_y)[0][i], 1e-14 );
        }
        else
        {
            EXPECT_SOFTEQ( 0.0, (*d_y)[0][i], 1e-14 );
        }
    }

    // New vector
    for( int i=0; i<d_y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        (*d_x)[0][i] = static_cast<double>(global_row+1);
    }

    d_operator->Apply(*d_x,*d_y);
    for( int i=0; i<d_y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        if( global_row == 19 )
        {
            EXPECT_SOFTEQ( 21.0, (*d_y)[0][i], 1e-14 );
        }
        else
        {
            EXPECT_SOFTEQ( 0.0, (*d_y)[0][i], 1e-14 );
        }
    }

    // Now set a shift
    d_operator->set_shift(0.5);
    d_x->PutScalar(1.0);
    d_operator->Apply(*d_x,*d_y);

    // Matlab computed reference
    double ref[] = {
        0.5000,  -1.0000,  -1.5000,  -2.0000,  -2.5000,  -3.0000,  -3.5000,
       -4.0000,  -4.5000,  -5.0000,  -5.5000,  -6.0000,  -6.5000,  -7.0000,
       -7.5000,  -8.0000,  -8.5000,  -9.0000,  -9.5000,  -9.0000 };

    for( int i=0; i<d_y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref[global_row], (*d_y)[0][i], 1e-14 );
    }


    // Different vector
    for( int i=0; i<d_y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        (*d_x)[0][i] = static_cast<double>(global_row+1);
    }
    d_operator->Apply(*d_x,*d_y);


    // Matlab computed reference
    double ref2[] = {
        -0.5000,   -2.0000,   -4.5000,   -8.0000,  -12.5000, -18.0000,
       -24.5000,  -32.0000,  -40.5000,  -50.0000,  -60.5000, -72.0000,
       -84.5000,  -98.0000, -112.5000, -128.0000, -144.5000, -162.0000,
      -180.5000, -179.0000 };
    for( int i=0; i<d_y->MyLength(); ++i )
    {
        int global_row = d_A->GRID(i);
        EXPECT_SOFTEQ( ref2[global_row], (*d_y)[0][i], 1e-14 );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstShiftedOperator.cc
//---------------------------------------------------------------------------//
