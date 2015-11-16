// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Device_Vector_device.t.hh
 * \author Seth R Johnson
 * \date   Thu Aug 01 11:33:12 2013
 * \brief  Device_Vector template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Device_Vector_device_t_hh
#define cuda_utils_Device_Vector_device_t_hh

#include "Device_Vector.hh"

#include <cuda_runtime.h>
#include <utility>

#include "Utils/harness/DBC.hh"
#include "Utils/comm/Logger.hh"
#include "Utils/utils/View_Field.hh"
#include "CudaDBC.hh"
#include "Host_Vector.hh"
#include "Stream.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with number of elements
 *
 * This does *not* perform any initialization, which might require a kernel
 * call.
 */
template<typename T>
Device_Vector<arch::Device,T>::Device_Vector(size_t count)
    : d_size(count)
    , d_data(NULL)
    , d_is_initialized(false)
{
    Require(count > 0);

    // Allocate memory
    this->allocate();

    Ensure(d_size == count);
    Ensure(d_data != NULL);
    Ensure(!d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Device_Vector<arch::Device,T>::Device_Vector(const_View_Field_t hostvec)
    : d_size(hostvec.size())
    , d_data(NULL)
{
    Require(hostvec.size() > 0);

    // Allocate memory
    this->allocate();

    // Copy memory to device
    this->assign(hostvec);

    Ensure(d_size == hostvec.size());
    Ensure(d_data != NULL);
    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Device_Vector<arch::Device,T>::Device_Vector(const Host_Vector_t& hostvec)
    : d_size(hostvec.size())
    , d_data(NULL)
{
    Require(hostvec.size() > 0);

    // Allocate memory
    this->allocate();

    // Copy memory to device
    this->assign(hostvec);

    Ensure(d_size == hostvec.size());
    Ensure(d_data != NULL);
    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a copy of data in another GPU vector.
 *
 * This is an \c explicit constructor because we don't want to be unwittingly
 * doing memory copy operations on the GPU.
 *
 * \note If DBC is on, an assertion will be thrown if the GPU vector being
 * copied is uninitialized.
 */
template<typename T>
Device_Vector<arch::Device,T>::Device_Vector(const This& rhs)
    : d_size(rhs.d_size)
    , d_data(NULL)
    , d_is_initialized(true)
{
    Require(rhs.is_initialized());

    // Allocate memory
    this->allocate();

    // Copy memory from the other device vector
    CudaCall(cudaMemcpy(d_data, rhs.d_data, d_size * sizeof(T),
                cudaMemcpyDeviceToDevice));

    Ensure(d_size == rhs.size());
    Ensure(d_data != NULL);
    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Allocate device memory and gracefully fail if necessary
 *
 */
template<typename T>
void Device_Vector<arch::Device,T>::allocate()
{
    Require(d_data == NULL);
    Require(d_size > 0);

    try
    {
        CudaCall(cudaMalloc(&d_data, d_size * sizeof(T)));
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

    Ensure(d_data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Free GPU memory on destruct
 */
template<typename T>
Device_Vector<arch::Device,T>::~Device_Vector()
{
    try
    {
        CudaCall(cudaFree(d_data));
    }
    catch (const profugus::assertion& e)
    {
        log(profugus::WARNING)
            << "Error: failed to free device data "
            << "at " << d_data << ": " << e.what();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory
 */
template<typename T>
void Device_Vector<arch::Device,T>::assign(const_View_Field_t hostvec)
{
    Require(hostvec.size() == size());
    Require(d_data);

    // Copy from host to device
    CudaCall(cudaMemcpy(d_data, &hostvec[0], hostvec.size() * sizeof(T),
                cudaMemcpyHostToDevice));

    // Set our state to initialized
    d_is_initialized = true;

    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory from a Host_Vector
 *
 * If the host vector uses "mapped" memory, then it's directly accessible on
 * the device, so we can use a device-to-device transfer (i.e. skipping the
 * front side bus entirely).
 */
template<typename T>
void Device_Vector<arch::Device,T>::assign(const Host_Vector_t& hostvec)
{
    Require(hostvec.size() == size());
    Require(d_data);

    if (hostvec.is_mapped())
    {
        // Copy from device-mapped pointer to device
        CudaCall(cudaMemcpy(d_data, hostvec.data(), hostvec.size() * sizeof(T),
                    cudaMemcpyDeviceToDevice));
    }
    else
    {
        // Copy from pinned host memory to device
        // We have to use the cpu_data accessor in case of WRITE_COMBINED
        // memory (which would prevent us from using brackets)
        CudaCall(cudaMemcpy(d_data, hostvec.cpu_data(),
                    hostvec.size() * sizeof(T),
                    cudaMemcpyHostToDevice));
    }

    // Set our state to initialized
    d_is_initialized = true;

    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory
 */
template<typename T>
void Device_Vector<arch::Device,T>::assign(const This& rhs)
{
    Require(rhs.size() == size());
    Require(d_data);

    // Copy from host to device
    CudaCall(cudaMemcpy(d_data, rhs.d_data, rhs.size() * sizeof(T),
                cudaMemcpyDeviceToDevice));

    // Set our state to initialized
    d_is_initialized = true;

    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory from a Host_Vector
 *
 * If the host vector uses "mapped" memory, then it's directly accessible on
 * the device, so we can use a device-to-device transfer (i.e. skipping the
 * front side bus entirely).
 */
template<typename T>
void Device_Vector<arch::Device,T>::assign_async(
        const Host_Vector_t& hostvec,
        Stream_t& stream)
{
    Require(hostvec.size() == size());
    Require(d_data);

    if (hostvec.is_mapped())
    {
        // Copy from device-mapped pointer to device
        CudaCall(cudaMemcpyAsync(d_data, hostvec.data(),
                    hostvec.size() * sizeof(T),
                    cudaMemcpyDeviceToDevice,
                    stream.handle()));
    }
    else
    {
        // Copy from pinned host memory to device
        CudaCall(cudaMemcpyAsync(d_data, hostvec.cpu_data(),
                    hostvec.size() * sizeof(T),
                    cudaMemcpyHostToDevice,
                    stream.handle()));
    }

    // Set our state to initialized
    d_is_initialized = true;

    Ensure(d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two same-sized device vectors
 */
template<typename T>
void Device_Vector<arch::Device,T>::swap(This& rhs)
{
    std::swap(d_data, rhs.d_data);
    std::swap(d_is_initialized, rhs.d_is_initialized);
}

//---------------------------------------------------------------------------//
template <typename T>
void Device_Vector<arch::Device,T>::to_host(profugus::View_Field<T> out) const
{
    Require(size() == out.size());
    Require(is_initialized());

    CudaCall(cudaMemcpy(&out[0], data(), size() * sizeof(T),
                cudaMemcpyDeviceToHost));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a device-to-host transfer for pinned memory
 *
 * If this is mapped, it will be faster (as a device-to-device transfer).
 */
template <typename T>
void Device_Vector<arch::Device,T>::to_host(Host_Vector<T>& out) const
{
    Require(size() == out.size());
    Require(is_initialized());

    if (out.is_mapped())
    {
        // Copy to device-mapped writable pointer
        CudaCall(cudaMemcpy(out.data(), data(), size() * sizeof(T),
                    cudaMemcpyDeviceToDevice));
    }
    else
    {
        CudaCall(cudaMemcpy(&out[0], data(), size() * sizeof(T),
                    cudaMemcpyDeviceToHost));
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a device-to-host asynchronous transfer
 */
template <typename T>
void Device_Vector<arch::Device,T>::to_host_async(
        Host_Vector<T>& out,
        Stream_t&       stream) const
{
    Require(size() == out.size());
    Require(is_initialized());

    if (out.is_mapped())
    {
        // Copy to device-mapped writable pointer
        CudaCall(cudaMemcpyAsync(out.data(), data(), size() * sizeof(T),
                    cudaMemcpyDeviceToDevice,
                    stream.handle()));
    }
    else
    {
        CudaCall(cudaMemcpyAsync(&out[0], data(), size() * sizeof(T),
                    cudaMemcpyDeviceToHost,
                    stream.handle()));
    }
}

//===========================================================================//

} // end namespace cuda

#endif // cuda_utils_Device_Vector_device_t_hh

//---------------------------------------------------------------------------//
//                   end of cuda_utils/Device_Vector_device.t.hh
//---------------------------------------------------------------------------//
