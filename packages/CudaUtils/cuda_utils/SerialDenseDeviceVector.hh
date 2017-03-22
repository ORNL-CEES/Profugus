//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/SerialDenseDeviceVector.hh
 * \author Stuart Slattery
 * \date   Tue Jan 5 2016
 * \brief  SerialDenseDeviceVector class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_SerialDenseDeviceVector_hh
#define cuda_utils_SerialDenseDeviceVector_hh

#include <Teuchos_Array.hpp>

#include "CudaMacros.hh"
#include "CudaDBC.hh"

namespace cuda_utils
{

//---------------------------------------------------------------------------//
/*!
 * \class SerialDenseDeviceVector
 * \brief Dense device vector for scattering data. On-device only. Does not
 * own memory.
 */
class SerialDenseDeviceVector
{
  public:

    // Device data constructor. Device-only for wrapping of data on-device.
    PROFUGUS_DEVICE_FUNCTION
    SerialDenseDeviceVector( const int size,
			     double* device_data )
	: d_size( size )
	, d_data( device_data )
    {
	PROFUGUS_INSIST_ON_DEVICE;
    }

    // Destructor. Prohibits copy construction and assignment. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    ~SerialDenseDeviceVector() { /* ... */ }

    // Get the number of rows. Host-accesible.
    PROFUGUS_DEVICE_FUNCTION
    int size() const
    {
	PROFUGUS_INSIST_ON_DEVICE;
	return d_size;
    }

    // Const value accessor. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    const double& operator()( const int i ) const
    {
	PROFUGUS_INSIST_ON_DEVICE;
	DEVICE_REQUIRE( i < d_size );
	return d_data[i];
    }

    // Non-const value accessor. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    double& operator()( const int i )
    {
	PROFUGUS_INSIST_ON_DEVICE;
	DEVICE_REQUIRE( i < d_size );
	return d_data[i];
    }

  private:

    // Length of vector.
    int d_size;

    // Vector data on the device.
    double* d_data;
};

//---------------------------------------------------------------------------//

} // end namespace cuda_utils

#endif // cuda_utils_SerialDenseDeviceVector_hh

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceVector.hh
//---------------------------------------------------------------------------//
