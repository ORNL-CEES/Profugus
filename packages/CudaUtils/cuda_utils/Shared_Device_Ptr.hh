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

#include "harness/DBC.hh"
#include "CudaDBC.hh"

namespace cuda
{
//===========================================================================//
/*!
 * \class Shared_Device_Ptr
 *
 * \brief Manage a pointer to an object on the device with a shared_ptr.
 *
 * Shared_Device_Ptr will manage a SHALLOW COPY of a host object on the
 * device. This means a host copy can manage device data automatically and we
 * can pass the device copy around to kernels without deallocating that
 * memory. When the device copy is done being used (i.e. this class destructor
 * is called) both the device copy will automatically be destroyed and the
 * host copy reference count will be decremented (or destroyed if zero). This
 * paradigm is most useful for managing C-style structs passed to CUDA kernels
 * (i.e. those with plain-old-data or raw pointers to device data).
 *
 * Shared_Device_Ptr may only be constructed from CUDA code as it is meant to
 * copy data to the device and manage the deallocation of that data on the
 * device.
 */
template<class T>
class Shared_Device_Ptr
{
  public:

    //! Default constructor.
    Shared_Device_Ptr()
    { /* ... */ }

    //! Host pointer constructor.
    Shared_Device_Ptr( const std::shared_ptr<T>& host_ptr )
	: d_host_ptr( host_ptr )
    {
	shallow_copy_host_object_to_device();
    }

    //! Get the shared pointer to host data managed by this object.
    inline std::shared_ptr<T>& get_host_ptr() { return d_host_ptr; }
    inline const std::shared_ptr<T>& get_host_ptr() const { return d_host_ptr; }

    //! Get the raw pointer to device data managed by this object. Accessible
    //! from the host or the device.
    T* get_device_ptr() { return d_device_ptr.get(); }
    inline const T* get_device_ptr() const { return d_device_ptr.get(); }

    //! Update device-side object from host object
    void update_device()
    {
#ifdef __NVCC__
        REQUIRE( d_host_ptr );
        REQUIRE( d_device_ptr );
        cudaMemcpy( d_device_ptr.get(), d_host_ptr.get(), sizeof(T),
                    cudaMemcpyHostToDevice );

        // No need to synchronize, no kernel launch on this stream can
        // start until transfer has completed
#endif
    }

    //! Update host-side object from device object
    void update_host()
    {
#ifdef __NVCC__
        REQUIRE( d_host_ptr );
        REQUIRE( d_device_ptr );
        cudaMemcpy( d_host_ptr.get(), d_device_ptr.get(), sizeof(T),
                    cudaMemcpyDeviceToHost );

        // Need to synchronize, transfer is asynchronous with respect to host
        cudaDeviceSynchronize();
#endif
    }


  private:

    // Smart pointer to HOST object.
    std::shared_ptr<T> d_host_ptr;

    // Smart pointer to the shallow copy of the object on the DEVICE.
    std::shared_ptr<T> d_device_ptr;

  private:

    // Shallow copy the host object to the device.
    void shallow_copy_host_object_to_device()
    {
#ifdef __NVCC__
	T* device_ptr;
	cudaMalloc( (void**) &device_ptr, sizeof(T) );
	cudaMemcpy( device_ptr, d_host_ptr.get(), sizeof(T),
		    cudaMemcpyHostToDevice );
	d_device_ptr =
	    std::shared_ptr<T>( device_ptr, [](T* t){ cudaFree(t); } );
#else
	DEVICE_INSIST(false,"Shared_Device_Ptr can only be constructed with NVCC!");
#endif // end __NVCC__
    }
};

//---------------------------------------------------------------------------//
// Free function for creating a Shared_Device_Ptr from object construction
// arguments.
template<class T,class ...Args>
Shared_Device_Ptr<T> shared_device_ptr( Args&&... args )
{
    return Shared_Device_Ptr<T>(
	std::make_shared<T>(std::forward<Args>(args)...) );
}

//===========================================================================//
} // end namespace cuda

#endif // cuda_utils_Shared_Device_Ptr_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Shared_Device_Ptr.hh
//---------------------------------------------------------------------------//
