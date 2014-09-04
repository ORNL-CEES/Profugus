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
#include "AnasaziOperatorTraits.hpp"

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
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

    // Initialization that are performed for each test
    void SetUp()
    {
        // Build an map
        d_N = 20;

        // Build CrsMatrix
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",d_N);
        d_B = linalg_traits::build_matrix<Matrix>("diagonal",d_N);

        // Build eigenvector
        d_x = linalg_traits::build_vector<MV>(d_N);
        d_y = linalg_traits::build_vector<MV>(d_N);

        // Build solver
        d_operator = Teuchos::rcp(new profugus::ShiftedOperator<MV,OP>());
        CHECK(!d_operator.is_null());
        d_operator->set_operator(d_A);
        d_operator->set_rhs_operator(d_B);
    }

  protected:

    int d_N;

    Teuchos::RCP<Matrix> d_A;
    Teuchos::RCP<Matrix> d_B;
    Teuchos::RCP<MV>     d_x;
    Teuchos::RCP<MV>     d_y;

    Teuchos::RCP<profugus::ShiftedOperator<MV,OP> > d_operator;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ShiftedOperatorTest, basic)
{
    // Unshifted operator, same as multiplying by A
    d_operator->set_shift(0.0);
    std::vector<double> one(d_N,1.0);
    linalg_traits::fill_vector<MV>(d_x,one);
    OPT::Apply(*d_operator,*d_x,*d_y);

    std::vector<double> ref(d_N,0.0);
    ref[0]  = 1.0;
    ref[19] = 1.0;
    linalg_traits::test_vector<MV>(d_y,ref);

    // New vector
    std::vector<double> vals(d_N);
    for( int i=0; i<d_N; ++i )
        vals[i] = static_cast<double>(i+1);
    linalg_traits::fill_vector<MV>(d_x,vals);

    OPT::Apply(*d_operator,*d_x,*d_y);
    std::fill(ref.begin(),ref.end(),0.0);
    ref[19] = 21.0;
    linalg_traits::test_vector<MV>(d_y,ref);

    // Now set a shift
    d_operator->set_shift(0.5);
    linalg_traits::fill_vector<MV>(d_x,one);
    OPT::Apply(*d_operator,*d_x,*d_y);

    // Matlab computed reference
    std::vector<double> ref2 = {
        0.5000,  -1.0000,  -1.5000,  -2.0000,  -2.5000,  -3.0000,  -3.5000,
       -4.0000,  -4.5000,  -5.0000,  -5.5000,  -6.0000,  -6.5000,  -7.0000,
       -7.5000,  -8.0000,  -8.5000,  -9.0000,  -9.5000,  -9.0000 };

    linalg_traits::test_vector<MV>(d_y,ref2);

    // Different vector
    linalg_traits::fill_vector<MV>(d_x,vals);
    OPT::Apply(*d_operator,*d_x,*d_y);

    // Matlab computed reference
    std::vector<double> ref3 = {
        -0.5000,   -2.0000,   -4.5000,   -8.0000,  -12.5000, -18.0000,
       -24.5000,  -32.0000,  -40.5000,  -50.0000,  -60.5000, -72.0000,
       -84.5000,  -98.0000, -112.5000, -128.0000, -144.5000, -162.0000,
      -180.5000, -179.0000 };

    linalg_traits::test_vector<MV>(d_y,ref3);
}

//---------------------------------------------------------------------------//
//                        end of tstShiftedOperator.cc
//---------------------------------------------------------------------------//
