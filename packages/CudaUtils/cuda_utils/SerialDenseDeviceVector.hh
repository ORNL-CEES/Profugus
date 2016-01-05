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

#include "harness/DBC.hh"
#include "CudaMacros.hh"

namespace cuda
{

//---------------------------------------------------------------------------//
/*!
 * \class SerialDenseDeviceVector
 * \brief Device matrix in ROW-MAJOR order for scattering data.
 */
class SerialDenseDeviceVector
{
  public:
    
    // Size constructor.
    SerialDenseDeviceVector( const int size, const double fill_value = 0.0 );

    // Host-data constructor.
    SerialDenseDeviceVector( const Teuchos::Array<double>& host_data );

    // Destructor. Prohibits copy construction and assignment.
    ~SerialDenseDeviceVector();

    // Get the number of rows. Host-accesible.
    PROFUGUS_HOST_DEVICE_FUNCTION(
	int size() const { return d_size; }
	)


    // Const value accessor. Device-only.
    PROFUGUS_DEVICE_FUNCTION(
	const double& operator()( const int i ) const
	{
	    return d_data[i];
	}
	)

    // Non-const value accessor. Device-only.
    PROFUGUS_DEVICE_FUNCTION(
	double& operator()( const int i )
	{
	    return d_data[i];
	}
	)
    
  private:

    // Allocate and copy host data to device.
    void allocate_and_copy( const Teuchos::Array<double>& host_data );

  private:

    // Length of vector.
    int d_size;

    // Vector data on the device.
    double* d_data;
};

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_SerialDenseDeviceVector_hh

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceVector.hh
//---------------------------------------------------------------------------//
