//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/SerialDenseMatrixVector_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_SerialDenseMatrixVector_Tester_hh
#define cuda_utils_test_SerialDenseMatrixVector_Tester_hh

#include <Teuchos_Array.hpp>
#include <Teuchos_TwoDArray.hpp>

//---------------------------------------------------------------------------//
// Matrix-Vector product. A*x = y
//---------------------------------------------------------------------------//
class SerialDenseMatrixVectorProduct
{
  public:

    SerialDenseMatrixVectorProduct( const Teuchos::Array<double>& A,
				    const Teuchos::Array<double>& x );

    ~SerialDenseMatrixVectorProduct();

    void multiply_kernel_launch();

    Teuchos::Array<double> get_result() const;

  private:

    int d_N;
    double* d_A;
    double* d_x;
    double* d_y;
};

//---------------------------------------------------------------------------//

#endif // cuda_utils_test_SerialDenseMatrixVector_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/SerialDenseMatrixVector_Tester.cuh
//---------------------------------------------------------------------------//
