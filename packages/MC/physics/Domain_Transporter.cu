//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Domain_Transporter.cu
 * \author Thomas M. Evans
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cuda.h>

#include "Domain_Transporter.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// Do something.
__global__ void do_something( int* values )
{
    int id = threadIdx.x + blockIdx.x * blockSize.x;
    values[id] = id;
}

//---------------------------------------------------------------------------//
void Domain_Transporter::launch_cuda() const
{
    int num_thread = 64;
    int num_block = 64;
    int num_data = num_thread*num_block;

    int* device_data;
    cudaMalloc( device_data, num_data*sizeof(int) );

    int* host_data = (int*) calloc( num_data, sizeof(int) );

    do_something<<<num_block,num_thread>>>( device_data );

    cudaMemcpy( host_data, device_data, num_data*sizeof(int),
		cudaMemcpyDeviceToHost );

    for ( int i = 0; i < num_data; ++i )
	std::cout << i << " " << host_data[i] << std::endl;

    cudaFree( device_data );
    free( host_data );
}

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.cu
//---------------------------------------------------------------------------//
