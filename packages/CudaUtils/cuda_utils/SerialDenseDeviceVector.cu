//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/SerialDenseDeviceVector.cu
 * \author Stuart Slattery
 * \date   Tue Jan 5 2016
 * \brief  SerialDenseDeviceVector class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SerialDenseDeviceVector.hh"

#include <cuda_runtime.h>

namespace cuda
{
//---------------------------------------------------------------------------//
// Size constructor.
SerialDenseDeviceVector::SerialDenseDeviceVector( const int size,
						  const double fill_value )
    : d_size( size )
{
    Teuchos::Array<double> host_data( size, fill_value );
    allocate_and_copy( host_data );
}

//---------------------------------------------------------------------------//
// Host-data constructor.
SerialDenseDeviceVector::SerialDenseDeviceVector( 
    const Teuchos::Array<double>& host_data )
    : d_size( host_data.size() )
{
    allocate_and_copy( host_data );
}

//---------------------------------------------------------------------------//
// Destructor. Prohibits copy construction and assignment.
SerialDenseDeviceVector::~SerialDenseDeviceVector()
{
    cudaFree( d_data );
}

//---------------------------------------------------------------------------//
// Allocate and copy host data to device.
void SerialDenseDeviceVector::allocate_and_copy( 
    const Teuchos::Array<double>& host_data )
{
    int mem_size = d_size*sizeof(double);
    cudaMalloc( (void**) &d_data, mem_size );
    cudaMemcpy( d_data, host_data.getRawPtr(), mem_size, 
		cudaMemcpyHostToDevice );
}

//---------------------------------------------------------------------------//

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceVector.cu
//---------------------------------------------------------------------------//
