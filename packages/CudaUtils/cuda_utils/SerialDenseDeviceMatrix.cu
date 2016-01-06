//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/SerialDenseDeviceMatrix.cu
 * \author Stuart Slattery
 * \date   Tue Jan 5 2016
 * \brief  SerialDenseDeviceMatrix class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SerialDenseDeviceMatrix.hh"

#include <cuda_runtime.h>

namespace cuda
{
//---------------------------------------------------------------------------//
// Size constructor.
SerialDenseDeviceMatrix::SerialDenseDeviceMatrix( const int num_rows,
						  const int num_cols,
						  const double fill_value )
    : d_num_rows( num_rows )
    , d_num_cols( num_cols )
    , d_owns_mem( true )
{
    Teuchos::TwoDArray<double> host_data( num_rows, num_cols, fill_value );
    allocate_and_copy( host_data );
}

//---------------------------------------------------------------------------//
// Host-data constructor.
SerialDenseDeviceMatrix::SerialDenseDeviceMatrix( 
    const int num_rows,
    const int num_cols,
    const Teuchos::Array<double>& host_data )
    : d_num_rows( num_rows )
    , d_num_cols( num_cols )
    , d_owns_mem( true )
{
    allocate_and_copy( host_data );
}

//---------------------------------------------------------------------------//
// Device data constructor.
SerialDenseDeviceMatrix::SerialDenseDeviceMatrix( const int num_rows, 
						  const int num_cols,
						  const bool owns_mem,
						  double* device_data )
    : d_num_rows( num_rows )
    , d_num_cols( num_cols )
    , d_owns_mem( owns_mem )
    , d_data( device_data )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor. Prohibits copy construction and assignment.
SerialDenseDeviceMatrix::~SerialDenseDeviceMatrix()
{
    if ( d_owns_mem )
    {
	cudaFree( d_data );
    }
}

//---------------------------------------------------------------------------//
// Allocate and copy host data to the device.
void SerialDenseDeviceMatrix::allocate_and_copy( 
    const Teuchos::TwoDArray<double>& host_data )
{
    int mem_size = d_num_rows*d_num_cols*sizeof(double);
    cudaMalloc( (void**) &d_data, mem_size );
    cudaMemcpy( d_data, host_data.getRawPtr(), mem_size,
		cudaMemcpyHostToDevice );
}

//---------------------------------------------------------------------------//

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceMatrix.cu
//---------------------------------------------------------------------------//
