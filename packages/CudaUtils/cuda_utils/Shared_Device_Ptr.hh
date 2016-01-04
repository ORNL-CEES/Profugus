//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Shared_Device_Ptr
 * \author Stuart Slattery, Tom Evans, Steven Hamilton
 * \date   Thu Dec 17 11:43:04 2015
 * \brief  Shared pointer to device data.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Shared_Device_Ptr_hh
#define cuda_utils_Shared_Device_Ptr_hh

#include <memory>

#include <cuda_runtime.h>

#include <harness/DBC.hh>

namespace cuda
{
//===========================================================================//
/*!
 * \class Shared_Device_Ptr
 *
 * \brief Manage a pointer to an object on the device with a shared_ptr.
 *
 * Shared_Device_Ptr will manage a COPY of a device object. It will not mirror
 * an equivalent host object.
 */
template<class T>
class Shared_Device_Ptr
{
  public:

    //! Default constructor.
    Shared_Device_Ptr()
    { /* ... */ }

    //! Copy constructor.
    Shared_Device_Ptr( const T& host_object )
    { 
	copy_host_object_to_device( host_object ); 
    }

    //! Arugment constructor. Arguments must be for a valid constructor of
    //! type T.
    template<class ...Args>
    explicit Shared_Device_Ptr( Args&&... args )
    {
	T host_object( std::forward<Args>(args)... );
	copy_host_object_to_device( host_object ); 
    }

    //! Get the raw pointer to device data managed by this object.
    inline T* get_device_ptr() { return d_device_ptr.get(); }
    inline const T* get_device_ptr() const { return d_device_ptr.get(); }

  private:

    // Smart pointer to DEVICE data.
    std::shared_ptr<T> d_device_ptr;

  private:

    // Copy a reference of T to the device.
    void copy_host_object_to_device( const T& host_object )
    {
#ifdef __NVCC__
	T* device_object;
	cudaMalloc( (void**) &device_object, sizeof(T) );
	cudaMemcpy( device_object, &host_object, sizeof(T),
		    cudaMemcpyHostToDevice );
	d_device_ptr = std::shared_ptr<T>( device_object, 
					   [](T* t){ CudaCall(cudaFree(t)); } );
#else
	INSIST( false, "Shared_Device_Ptr can only be constructed with NVCC!" );
#endif // end __NVCC__
    }
};

//===========================================================================//
} // end namespace cuda

#endif // cuda_utils_Shared_Device_Ptr_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Shared_Device_Ptr.hh
//---------------------------------------------------------------------------//
