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
__global__ get_vector_kernel( const cuda_profugus::XS_Device* xs,
			      const int matid,
			      const int type,
			      double* out_vector )
{
    cuda::SerialDenseDeviceVector in_vector = xs->vector( matid, type );
    int i = threadIdx.x;
    out_vector[i] = in_vector(i);
}

//---------------------------------------------------------------------------//
__global__ get_matrix_kernel( const cuda_profugus::XS_Device* xs,
			      const int matid,
			      const int pn,
			      double* out_matrix )
{
    cuda::SerialDenseDeviceMatrix in_matrix = xs->matrix( matid, pn );
    int i = threadIdx.x;
    out_matrix[i] = in_matrix(i);
}

//---------------------------------------------------------------------------//
// XS_Device_Tester
//---------------------------------------------------------------------------//
XS_Device_Tester::XS_Device_Tester( const profugus::XS& xs )
    : d_xs( xs )
{ /* ... */ }

//---------------------------------------------------------------------------//
const typename profugus::XS::Vector& 
XS_Device_Tester::vector( int matid, int type ) const
{
    int N = d_xs.get_host_ptr()->num_group();

    double* device_vec;
    cudaMalloc( (void**) &device_vec, N * sizeof(double) );

    get_vector_kernel<<<1,N>>>( d_xs.get_device_ptr(),
				matid,
				type,
				device_vec );

    profugus::XS::Vector host_vec( N );
    cudaMemcpy( host_vec.values, device_vec.values(), N * sizeof(double),
		cudaMemcpyDeviceToHost );
    cudaFree( device_vec );
    return host_vec;
}

//---------------------------------------------------------------------------//
const typename profugus::XS::Matrix& 
XS_Device_Tester::matrix( int matid, int pn ) const
{
    int N = d_xs.get_host_ptr()->num_group() *
	    d_xs.get_host_ptr()->num_group();

    double* device_mat;
    cudaMalloc( (void**) &device_mat, N * sizeof(double) );

    get_matrix_kernel<<<1,N>>>( d_xs.get_device_ptr(),
				matid,
				type,
				device_mat );

    profugus::XS::Matrix host_mat( N, N );
    cudaMemcpy( host_mat.values, device_mat.values(), N * sizeof(double),
		cudaMemcpyDeviceToHost );
    cudaFree( device_mat );
    return host_mat;
}

//---------------------------------------------------------------------------//
//                 end of cuda_utils/XS_Device_Tester.cu
//---------------------------------------------------------------------------//
