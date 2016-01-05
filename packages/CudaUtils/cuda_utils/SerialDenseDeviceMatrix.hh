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

#include <Teuchos_TwoDArray.hpp>

#include "harness/DBC.hh"

namespace cuda
{

//---------------------------------------------------------------------------//
/*!
 * \class SerialDenseDeviceMatrix
 * \brief Device matrix in ROW-MAJOR order for scattering data.
 */
class SerialDenseDeviceMatrix
{
  public:
    
    // Size constructor.
    SerialDenseDeviceMatrix( const int num_rows, const int num_cols,
			     const double fill_value = 0.0 );

    // Host data constructor.
    SerialDenseDeviceMatrix( const Teuchos::TwoDArray<double>& host_data );

    // Destructor. Prohibits copy construction and assignment.
    ~SerialDenseDeviceMatrix();

    // Get the number of rows. Host-accesible.
#ifdef __NVCC__
    __host__ __device__ 
#endif
    int num_rows() const { return d_num_rows; }

    // Get the number of columns. Host-accessible.
#ifdef __NVCC__
    __host__ __device__ 
#endif
    int num_cols() const { return d_num_cols; }

    // Const value accessor. Device-only.
#ifdef __NVCC__
    __device__ 
#endif
    const double& operator()( const int row, const int col ) const
    {
#ifdef __NVCC__
	return d_data[row*d_num_cols + col];
#else
	INSIST( false, "Tried to access matrix data on host!" );
	return d_data[0];
#endif
    }

    // Non-const value accessor. Device-only.
#ifdef __NVCC__
    __device__ 
#endif
    double& operator()( const int row, const int col )
    {
#ifdef __NVCC__
	return d_data[row*d_num_cols + col];
#else
	INSIST( false, "Tried to access matrix data on host!" );
	return d_data[0];
#endif
    }
    
  private:

    // Allocate and copy host data to the device.
    void allocate_and_copy( const Teuchos::TwoDArray<double>& host_data );

  private:

    // Number of rows.
    int d_num_rows;
    
    // Number of columns.
    int d_num_cols;

    // Matrix data.
    double* d_data;
};

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_SerialDenseDeviceMatrix_hh

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceMatrix.hh
//---------------------------------------------------------------------------//
