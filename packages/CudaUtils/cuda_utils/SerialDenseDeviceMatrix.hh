//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/SerialDenseDeviceMatrix.hh
 * \author Stuart Slattery
 * \date   Tue Jan 5 2016
 * \brief  SerialDenseDeviceMatrix class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_SerialDenseDeviceMatrix_hh
#define cuda_utils_SerialDenseDeviceMatrix_hh

#include <Teuchos_Array.hpp>

#include "CudaDBC.hh"
#include "CudaMacros.hh"

namespace cuda_utils
{

//---------------------------------------------------------------------------//
/*!
 * \class SerialDenseDeviceMatrix
 * \brief Device matrix in COLUMN-MAJOR order for scattering data. This
 * manages a view of the data and does not own the data.
 */
class SerialDenseDeviceMatrix
{
  public:
    
    // Device data constructor. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    SerialDenseDeviceMatrix( const int num_rows, 
			     const int num_cols,
			     double* device_data )
	: d_num_rows( num_rows )
	, d_num_cols( num_cols )
	, d_data( device_data )
    { 
	PROFUGUS_INSIST_ON_DEVICE;
    }

    // Destructor. Prohibits copy construction and assignment. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    ~SerialDenseDeviceMatrix() { /* ... */ }

    // Get the number of rows. Host-accesible.
    PROFUGUS_DEVICE_FUNCTION
    int num_rows() const 
    { 
	PROFUGUS_INSIST_ON_DEVICE;
	return d_num_rows; 
    }

    // Get the number of columns. Host-accessible.
    PROFUGUS_DEVICE_FUNCTION
    int num_cols() const 
    { 
	PROFUGUS_INSIST_ON_DEVICE;
	return d_num_cols; 
    }

    // Const value accessor. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    const double& operator()( const int row, const int col ) const
    {
	PROFUGUS_INSIST_ON_DEVICE;
	REQUIRE( row < d_num_rows );
	REQUIRE( col < d_num_cols );
	return d_data[col*d_num_rows + row];
    }

    // Non-const value accessor. Device-only.
    PROFUGUS_DEVICE_FUNCTION
    double& operator()( const int row, const int col )
    {
	PROFUGUS_INSIST_ON_DEVICE;
	REQUIRE( row < d_num_rows );
	REQUIRE( col < d_num_cols );
	return d_data[col*d_num_rows + row];
    }
    
  private:

    // Number of rows.
    int d_num_rows;
    
    // Number of columns.
    int d_num_cols;

    // Matrix data on the device.
    double* d_data;
};

//---------------------------------------------------------------------------//

} // end namespace cuda_utils

#endif // cuda_utils_SerialDenseDeviceMatrix_hh

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceMatrix.hh
//---------------------------------------------------------------------------//
