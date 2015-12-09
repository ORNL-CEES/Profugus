// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Hardware.cc
 * \author Seth R Johnson
 * \date   Tue Jul 09 15:43:18 2013
 * \brief  Device member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Hardware.hh"

#include <config.h>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#include "CudaDBC.hh"
#endif

#include <iostream>
#include "comm/global.hh"

using cuda::arch::Device;

namespace cuda
{
//===========================================================================//
// DEVICE SPECIALIZATION
//---------------------------------------------------------------------------//
// STATIC VARIABLES
//---------------------------------------------------------------------------//

//! Device ID after Acquire; -1 if unset.
int Hardware<Device>::d_device_id = -1;

//! Device name after Acquire
std::string Hardware<Device>::d_name = "n/a";

//! Number of SMs in the device
unsigned int Hardware<Device>::d_multiProcessorCount = 0;

//! Amount of constant memory available
unsigned int Hardware<Device>::d_totalConstMem = 0;

//! Amount of shared memory available in each SM (divisible among blocks)
unsigned int Hardware<Device>::d_sharedMemPerBlock = 0;

//! Number of cores per multiprocessor in the device
unsigned int Hardware<Device>::d_cores_per_mp = -1;

//! Number of cores per multiprocessor in the device
unsigned int Hardware<Device>::d_maxThreadsPerMultiProcessor = 0;

//! Default block size for kernel launches
unsigned int Hardware<Device>::d_default_block_size = 256;

//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Return whether a valid device exists.
 *
 * If CUDA support is not enabled, this will always return false
 */
bool Hardware<Device>::valid_device_exists()
{
    int num_devices = -1;
#ifdef USE_CUDA
    CudaCall(cudaGetDeviceCount(&num_devices));
#endif
    return num_devices > 0;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Choose a CUDA device based on some heuristic.
 */
void Hardware<Device>::acquire(Acquire_Method method)
{
    INSIST(d_device_id == -1, "A device has already been acquired.");
    VALIDATE(method < END_ACQUIRE_METHOD, "Invalid method " << method);

#ifdef USE_CUDA
    cudaDeviceProp prop;
    int num_devices = -1;
    CudaCall(cudaGetDeviceCount(&num_devices));

    int best_value = -1;
    int this_value = -1;
    for (int device = 0; device < num_devices; ++device)
    {
        CudaCall(cudaGetDeviceProperties(&prop, device));

        // Check whether this device is valid
        if (prop.computeMode == cudaComputeModeProhibited)
            continue;

        switch (method)
        {
            case HIGHEST_FLOPS:
                this_value = prop.multiProcessorCount
                             * num_cores_per_mp(prop.major, prop.minor)
                             * prop.clockRate;
                break;
            case HIGHEST_COMPUTE:
                this_value = (prop.major << 4) + prop.minor;
                break;
            default: // squelch warning
                break;
        }

        if (this_value > best_value)
        {
            best_value  = this_value;
            // Set device ID and store properties
            d_device_id           = device;
            d_name                = prop.name;
            d_multiProcessorCount = prop.multiProcessorCount;
            d_totalConstMem       = prop.totalConstMem;
            d_sharedMemPerBlock   = prop.sharedMemPerBlock;
            d_cores_per_mp        = num_cores_per_mp(prop.major, prop.minor);

            d_maxThreadsPerMultiProcessor = prop.maxThreadsPerMultiProcessor;
        }
    }

    INSIST(d_device_id != 1, "No valid devices were found.");
    CudaCall(cudaSetDevice(d_device_id));

    // Reset the device
    CudaCall(cudaDeviceReset());

    // Initialize device, causing the ~.2 second penalty to appear in this
    // function rather than in later initialization
    CudaCall(cudaFree(0));

    if (profugus::node() == 0)
    {
        std::size_t free, total;
        CudaCall(cudaMemGetInfo(&free, &total));

        std::cout << "*** Initialized device '"
            << d_name
            << "' (Device ID " << d_device_id << ")"
            << " which has " << (free >> 20) << " of "
            << (total >> 20) << "MB free."
            << std::endl;
    }


#else // USE_CUDA disabled
    INSIST(false, "CUDA is disabled in this build.");
#endif // USE_CUDA
}

//---------------------------------------------------------------------------//
/*!
 * \brief Ensure asynchronous copies/kernel launches are complete
 *
 * This will raise a profugus::assertion if any Cuda issues are encountered.
 */
void Hardware<Device>::synchronize()
{
#ifdef USE_CUDA
    CudaCall(cudaDeviceSynchronize());
#else // USE_CUDA disabled
    INSIST(false, "CUDA is disabled in this build.");
#endif //USE_CUDA
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of cores per multiprocessor
 *
 * This is derived from information in helper_cuda.h in the cuda samples.
 */
unsigned int Hardware<Device>::num_cores_per_mp(int major, int minor)
{
    const int val = (major << 4) + minor;
    switch (val)
    {
        case 0x10: return   8; // Tesla  Generation (SM 1.0) G80   class
        case 0x11: return   8; // Tesla  Generation (SM 1.1) G8x   class
        case 0x12: return   8; // Tesla  Generation (SM 1.2) G9x   class
        case 0x13: return   8; // Tesla  Generation (SM 1.3) GT200 class
        case 0x20: return  32; // Fermi  Generation (SM 2.0) GF100 class
        case 0x21: return  48; // Fermi  Generation (SM 2.1) GF10x class
        case 0x30: return 192; // Kepler Generation (SM 3.0) GK10x class
        case 0x35: return 192; // Kepler Generation (SM 3.5) GK11x class
    }
    // Unknown SM version
    return 192;
}

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of Hardware.cc
//---------------------------------------------------------------------------//
