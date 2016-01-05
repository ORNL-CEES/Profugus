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

#include "../cuda_utils/Shared_Device_Ptr.hh"
#include "../cuda_utils/SerialDenseDeviceMatrix.hh"
#include "../cuda_utils/SerialDenseDeviceVector.hh"

//---------------------------------------------------------------------------//
// Matrix-Vector product. A*x = y
//---------------------------------------------------------------------------//
class SerialDenseMatrixVectorProduct
{
  public:

    SerialDenseMatrixVectorProduct( const Teuchos::TwoDArray<double>& A,
				    const Teuchos::Array<double>& x );

    void multiply_kernel_launch();

    Teuchos::Array<double> get_result() const;

  private:

    cuda::Shared_Device_Ptr<cuda::SerialDenseDeviceMatrix> d_A;
    cuda::Shared_Device_Ptr<cuda::SerialDenseDeviceVector> d_x;
    cuda::Shared_Device_Ptr<cuda::SerialDenseDeviceVector> d_y;
};

//---------------------------------------------------------------------------//

#endif // cuda_utils_test_SerialDenseMatrixVector_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/SerialDenseMatrixVector_Tester.cuh
//---------------------------------------------------------------------------//
