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

#include <exception>

#include "CudaDBC.hh"
#include "Stream.hh"

#include "comm/Logger.hh"

#include <cuda_runtime.h>

namespace cuda_utils
{

namespace memory
{
//---------------------------------------------------------------------------//
// Memory Management Functions.
//---------------------------------------------------------------------------//
// Cuda allocate memory. Will allocate memory on the device of the requested
// size and check for errors.
template<typename T>
void Malloc( T*& device_ptr, const std::size_t N )
{
    try
    {
	CudaCall( cudaMalloc( (void**) &device_ptr, N * sizeof(T) ) );
    }
    catch (const std::exception& e)
    {
        try
        {
            // Give more info on memory needs and availability
            std::size_t free, total;
            CudaCall(cudaMemGetInfo(&free, &total));
            profugus::log(profugus::WARNING)
                << "Error: Failed to allocate "
                << N * sizeof(T) << " bytes on device: only "
                << free << " of " << total << " bytes are free";
        }
        catch (const std::exception& e)
        {
	    profugus::log(profugus::WARNING)
                << "Error: Failed to allocate on device "
                << "and failed to get information about the failure.";
        }

        throw e;
    }
}

//---------------------------------------------------------------------------//
// Cuda free memory. This may be called in a destructor so we check for errors
// but log them instead of re-throwing the exception.
template<typename T>
void Free( T* device_ptr )
{
    try
    {
	CudaCall( cudaFree(device_ptr) );
    }
    catch (const std::exception& e)
    {
	profugus::log(profugus::WARNING)
	    << "Error: failed to free device data "
	    << "at " << device_ptr << ": " << e.what();
    }
}

//---------------------------------------------------------------------------//
// Copy from the host to the device.
template<typename T>
void Copy_To_Device( T* device_ptr, 
		     const T* host_ptr, 
		     const std::size_t N )
{
    CudaCall( cudaMemcpy(device_ptr, host_ptr, N*sizeof(T),
			 cudaMemcpyHostToDevice) );
}

//---------------------------------------------------------------------------//
// Copy from the device to the host.
template<typename T>
void Copy_To_Host( T* host_ptr, 
		   const T* device_ptr, 
		   const std::size_t N )
{
    CudaCall( cudaMemcpy(host_ptr, device_ptr, N*sizeof(T),
			 cudaMemcpyDeviceToHost) );
}

//---------------------------------------------------------------------------//
// Asynchronous copy from the host to the device.
template<typename T>
void Copy_To_Device_Async( T* device_ptr, 
			   const T* host_ptr, 
			   const std::size_t N,
			   ::cuda_utils::Stream<::cuda_utils::arch::Device>& stream )
{
    CudaCall( cudaMemcpyAsync(device_ptr, host_ptr, N*sizeof(T),
			      cudaMemcpyHostToDevice, stream.handle()) );
}

//---------------------------------------------------------------------------//
// Asynchronous copy from the device to the host.
template<typename T>
void Copy_To_Host_Async( T* host_ptr, 
			 const T* device_ptr, 
			 const std::size_t N,
			 ::cuda_utils::Stream<::cuda_utils::arch::Device>& stream )
{
    CudaCall( cudaMemcpyAsync(host_ptr, device_ptr, N*sizeof(T),
			      cudaMemcpyDeviceToHost, stream.handle()) );
}

//---------------------------------------------------------------------------//

} // end namespace memory

} // end namespace cuda_utils

#endif // CudaUtils_cuda_utils_Memory_cuh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Utility_Functions.hh
//---------------------------------------------------------------------------//
