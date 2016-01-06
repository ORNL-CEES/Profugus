//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstSerialDenseMatrixVector.hh
 * \author Stuart Slattery
 * \date   Thu Dec 17 11:43:04 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include <Teuchos_TwoDArray.hpp>
#include <Teuchos_Array.hpp>

#include "SerialDenseMatrixVector_Tester.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class SerialDenseMatrixVectorTest : public ::testing::Test
{

  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(SerialDenseMatrixVectorTest, pointer_test)
{
    int N = 1000;

    // Make a matrix with COLUMN-MAJOR ORDER.
    Teuchos::Array<double> A( N*N );
    for ( int i = 0; i < N; ++i )
    {
	for ( int j = 0; j < N; ++j )
	{
	    A[j*N + i] = i;
	}
    }

    // Make a vector.
    Teuchos::Array<double> x( N );
    for ( int i = 0; i < N; ++i )
    {
	x[i] = i;
    }

    // Create a matrix-vector product tester.
    SerialDenseMatrixVectorProduct product( A, x );

    // Compute the matrix-vector product on the GPU.
    product.multiply_kernel_launch();

    // Calculate the expected result.
    double gold = 0.0;
    for ( int i = 0; i < N; ++i )
    {
	gold += x[i];
    }

    // Get the result on the host and check it.
    Teuchos::Array<double> result = product.get_result();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( i*gold, result[i] );
    }    
}

//---------------------------------------------------------------------------//
//                        end of tstSerialDenseMatrixVector.cc
//---------------------------------------------------------------------------//
