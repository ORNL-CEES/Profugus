//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/SerialDenseMatrixVector_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SerialDenseMatrixVector_Tester.hh"

#include "../cuda_utils/SerialDenseDeviceMatrix.hh"
#include "../cuda_utils/SerialDenseDeviceVector.hh"

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// Matrix-Vector product kernel.
__global__ void multiply_kernel( const int N,
				 double* A,
				 double* x,
				 double* y )
{
    const cuda_utils::SerialDenseDeviceMatrix Amat( N, N, A );
    const cuda_utils::SerialDenseDeviceVector xvec( N , x );
    cuda_utils::SerialDenseDeviceVector yvec( N, y );

    int i = threadIdx.x;
    yvec(i) = 0.0;
    for ( int j = 0; j < Amat.num_cols(); ++j )
    {
	yvec(i) += Amat(i,j) * xvec(j);
    }
}

//---------------------------------------------------------------------------//
// Product result copy kernel.
__global__ void result_copy_kernel( const double* y,
				    double* result )
{
    int i = threadIdx.x;
    result[i] = y[i];
}

//---------------------------------------------------------------------------//
// SerialDenseMatrixVectorProduct class functions.
SerialDenseMatrixVectorProduct::SerialDenseMatrixVectorProduct( 
    const Teuchos::Array<double>& A,
    const Teuchos::Array<double>& x )
{
    d_N = x.size();

    cudaMalloc( (void**) &d_A, A.size() * sizeof(double) );
    cudaMemcpy( d_A, A.getRawPtr(), A.size() * sizeof(double),
		cudaMemcpyHostToDevice );

    cudaMalloc( (void**) &d_x, x.size() * sizeof(double) );
    cudaMemcpy( d_x, x.getRawPtr(), x.size() * sizeof(double),
		cudaMemcpyHostToDevice );

    cudaMalloc( (void**) &d_y, x.size() * sizeof(double) );
}

//---------------------------------------------------------------------------//
SerialDenseMatrixVectorProduct::~SerialDenseMatrixVectorProduct()
{
    cudaFree( d_A );
    cudaFree( d_x );
    cudaFree( d_y );
}

//---------------------------------------------------------------------------//
void SerialDenseMatrixVectorProduct::multiply_kernel_launch()
{
    multiply_kernel<<<1,d_N>>>( d_N, d_A, d_x, d_y );
}

//---------------------------------------------------------------------------//
Teuchos::Array<double> SerialDenseMatrixVectorProduct::get_result() const
{
    int mem_size = d_N*sizeof(double);

    double* device_result;
    cudaMalloc( (void**) &device_result, mem_size );
    result_copy_kernel<<<1,d_N>>>( d_y, device_result );

    Teuchos::Array<double> host_result( d_N );
    cudaMemcpy( host_result.getRawPtr(), device_result, mem_size,
		cudaMemcpyDeviceToHost );

    cudaFree( device_result );

    return host_result;
}

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/SerialDenseMatrixVector_Tester.cuh
//---------------------------------------------------------------------------//
