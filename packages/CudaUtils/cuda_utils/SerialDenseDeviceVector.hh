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
#ifdef __NVCC__
    __host__ __device__ 
#endif
    int size() const { return d_size; }

    // Const value accessor. Device-only.
#ifdef __NVCC__
    __device__ 
#endif
    const double& operator()( const int i ) const
    {
#ifdef __NVCC__
	return d_data[i];
#else
	INSIST( false, "Tried to access vector data on host!" );
	return d_data[0];
#endif
    }

    // Non-const value accessor. Device-only.
#ifdef __NVCC__
    __device__ 
#endif
    double& operator()( const int i )
    {
#ifdef __NVCC__
	return d_data[i];
#else
	INSIST( false, "Tried to access vector data on host!" );
	return d_data[0];
#endif
    }
    
  private:

    // Allocate and copy host data to device.
    void allocate_and_copy( const Teuchos::Array<double>& host_data );

  private:

    // Length of vector.
    int d_size;

    // Vector data.
    double* d_data;
};

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_SerialDenseDeviceVector_hh

//---------------------------------------------------------------------------//
//                 end of SerialDenseDeviceVector.hh
//---------------------------------------------------------------------------//
