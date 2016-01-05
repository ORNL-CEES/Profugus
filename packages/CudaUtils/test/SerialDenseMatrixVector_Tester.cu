//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/SerialDenseMatrixVector_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SerialDenseMatrixVector_Tester.hh"

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// Matrix-Vector product kernel.
__global__ void multiply_kernel( const cuda::SerialDenseDeviceMatrix* A,
				 const cuda::SerialDenseDeviceVector* x,
				 cuda::SerialDenseDeviceVector* y )
{
    int i = threadIdx.x;
    (*y)(i) = 0.0;
    for ( int j = 0; j < A->num_cols(); ++j )
    {
	(*y)(i) += (*A)(i,j) * (*x)(j);
    }
}

//---------------------------------------------------------------------------//
// Product result copy kernel.
__global__ void result_copy_kernel( const cuda::SerialDenseDeviceVector* y,
				    double* result )
{
    int i = threadIdx.x;
    result[i] = (*y)(i);
}

//---------------------------------------------------------------------------//
// SerialDenseMatrixVectorProduct class functions.
SerialDenseMatrixVectorProduct::SerialDenseMatrixVectorProduct( 
    const Teuchos::TwoDArray<double>& A,
    const Teuchos::Array<double>& x )
    : d_A( A )
    , d_x( x )
    , d_y( x.size() )
{ /* ... */ }

//---------------------------------------------------------------------------//
void SerialDenseMatrixVectorProduct::multiply_kernel_launch()
{
    int N = d_x.get_host_ptr()->size();
    multiply_kernel<<<1,N>>>( d_A.get_device_ptr(),
			      d_x.get_device_ptr(),
			      d_y.get_device_ptr() );
}

//---------------------------------------------------------------------------//
Teuchos::Array<double> SerialDenseMatrixVectorProduct::get_result() const
{
    int N = d_y.get_host_ptr()->size();
    int mem_size = N*sizeof(double);

    double* device_result;
    cudaMalloc( (void**) &device_result, mem_size );
    result_copy_kernel<<<1,N>>>( d_y.get_device_ptr(), device_result );

    Teuchos::Array<double> host_result( N );
    cudaMemcpy( host_result.getRawPtr(), device_result, mem_size,
		cudaMemcpyDeviceToHost );

    cudaFree( device_result );

    return host_result;
}

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/SerialDenseMatrixVector_Tester.cuh
//---------------------------------------------------------------------------//
