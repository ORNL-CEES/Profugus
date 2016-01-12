//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Memory.cuh
 * \author Stuart Slattery
 * \brief  Free function helpers for memory management.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Memory_cuh
#define CudaUtils_cuda_utils_Memory_cuh

#include "CudaDBC.hh"
#include "Stream.hh"

#include <cuda_runtime.h>

namespace cuda
{

namespace memory
{
//---------------------------------------------------------------------------//
// Memory Management Functions.
//---------------------------------------------------------------------------//
// Cuda malloc memory.
template<typename T>
void Malloc( T* device_ptr, const std::size_t N )
{
    try
    {
	CudaCall( cudaMalloc( (void**) &device_ptr, N * sizeof(T) ) );
    }
    catch (const profugus::assertion& e)
    {
        try
        {
            // Give more info on memory needs and availability
            std::size_t free, total;
            CudaCall(cudaMemGetInfo(&free, &total));
            profugus::log(profugus::WARNING)
                << "Error: Failed to allocate "
                << d_size * sizeof(T) << " bytes on device: only "
                << free << " of " << total << " bytes are free";
        }
        catch (const profugus::assertion& e)
        {
            log(profugus::WARNING)
                << "Error: Failed to allocate on device "
                << "and failed to get information about the failure.";
        }

        throw e;
    }
}

//---------------------------------------------------------------------------//
// Cuda free memory. This may be called in a destructor.
void Free( T* device_ptr )
{
    try
    {
        CudaCall( cudaFree(device_ptr) );
    }
    catch (const profugus::assertion& e)
    {
        log(profugus::WARNING)
            << "Error: failed to free device data "
            << "at " << d_data << ": " << e.what();
    }
}

//---------------------------------------------------------------------------//
// Copy from the host to the device.
template<typename T>
void Memcpy_To_Device( T* device_ptr, 
		       const T* host_ptr, 
		       const std::size_t N )
{
    cudaCall( cudaMemcpy(device_ptr, host_ptr, N*sizeof(T),
			 cudaMemcpyToHostToDevice) );
}

//---------------------------------------------------------------------------//
// Copy from the device to the host.
template<typename T>
void Memcpy_To_Host( T* host_ptr, 
		     const T* device_ptr, 
		     const std::size_t N )
{
    cudaCall( cudaMemcpy(host_ptr, device_ptr, N*sizeof(T),
			 cudaMemcpyToDeviceToHost) );
}

//---------------------------------------------------------------------------//
// Asynchronous copy from the host to the device.
template<typename T>
void Memcpy_To_Device_Async( T* device_ptr, 
			     const T* host_ptr, 
			     const std::size_t N,
			     ::cuda::Stream<::cuda::arch::Device>& stream )
{
    cudaCall( cudaMemcpyAsync(device_ptr, host_ptr, N*sizeof(T),
			      cudaMemcpyToHostToDevice, stream.handle()) );
}

//---------------------------------------------------------------------------//
// Asynchronous copy from the device to the host.
template<typename T>
void Memcpy_To_Host_Async( T* host_ptr, 
			   const T* device_ptr, 
			   const std::size_t N,
			   ::cuda::Stream<::cuda::arch::Device>& stream )
{
    cudaCall( cudaMemcpyAsync(host_ptr, device_ptr, N*sizeof(T),
			      cudaMemcpyToDeviceToHost, stream.handle()) );
}

//---------------------------------------------------------------------------//

} // end namespace memory

} // end namespace cuda

#endif // CudaUtils_cuda_utils_Memory_cuh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Utility_Functions.hh
//---------------------------------------------------------------------------//
