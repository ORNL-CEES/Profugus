//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_xs/test/XS_Device_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "XS_Device_Tester.hh"

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// CUDA Kernels
//---------------------------------------------------------------------------//
__global__ void get_vector_kernel( const cuda_profugus::XS_Device* xs,
				   const int matid,
				   const int type,
				   double* out_vector )
{
    cuda_utils::SerialDenseDeviceVector in_vector = xs->vector( matid, type );
    int i = threadIdx.x;
    out_vector[i] = in_vector(i);
}

//---------------------------------------------------------------------------//
__global__ void get_matrix_kernel( const cuda_profugus::XS_Device* xs,
				   const int matid,
				   const int pn,
				   double* out_matrix )
{
    cuda_utils::SerialDenseDeviceMatrix in_matrix = xs->matrix( matid, pn );
    int row = threadIdx.x;
    int col = threadIdx.y;
    int num_rows = in_matrix.num_rows();
    out_matrix[col*num_rows + row] = in_matrix(row,col);
}

//---------------------------------------------------------------------------//
// XS_Device_Tester
//---------------------------------------------------------------------------//
XS_Device_Tester::XS_Device_Tester( const profugus::XS& xs )
: d_xs( cuda_utils::shared_device_ptr<cuda_profugus::XS_Device>(xs) )
{ /* ... */ }

//---------------------------------------------------------------------------//
const typename profugus::XS::Vector
XS_Device_Tester::vector( int matid, int type ) const
{
    int N = d_xs.get_host_ptr()->num_groups();

    double* device_vec;
    cudaMalloc( (void**) &device_vec, N * sizeof(double) );

    get_vector_kernel<<<1,N>>>( d_xs.get_device_ptr(),
				matid,
				type,
				device_vec );

    profugus::XS::Vector host_vec( N );
    cudaMemcpy( host_vec.values(), device_vec, N * sizeof(double),
		cudaMemcpyDeviceToHost );
    cudaFree( device_vec );
    return host_vec;
}

//---------------------------------------------------------------------------//
const typename profugus::XS::Matrix 
XS_Device_Tester::matrix( int matid, int pn ) const
{
    int N = d_xs.get_host_ptr()->num_groups();

    double* device_mat;
    cudaMalloc( (void**) &device_mat, N * N * sizeof(double) );

    dim3 grid_dim(N,N);
    get_matrix_kernel<<<1,grid_dim>>>( d_xs.get_device_ptr(),
				       matid,
				       pn,
				       device_mat );

    profugus::XS::Matrix host_mat( N, N );
    cudaMemcpy( host_mat.values(), device_mat, N * N * sizeof(double),
		cudaMemcpyDeviceToHost );
    cudaFree( device_mat );
    return host_mat;
}

//---------------------------------------------------------------------------//
//                 end of cuda_utils/XS_Device_Tester.cu
//---------------------------------------------------------------------------//
