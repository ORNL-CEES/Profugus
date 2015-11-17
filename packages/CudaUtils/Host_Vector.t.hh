// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Host_Vector.t.hh
 * \author Seth R Johnson
 * \date   Mon Aug 12 08:48:53 2013
 * \brief  Host_Vector template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Host_Vector_t_hh
#define cuda_utils_Host_Vector_t_hh

#include "Host_Vector.hh"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "config.h"
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#include "comm/Logger.hh"
#include "utils/View_Field.hh"
#include "CudaDBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
/*!
 * \brief Allocate and initialize the host vector
 */
template<class T>
Host_Vector<T>::Host_Vector(
            size_t count,
            const_reference value,
            Alloc_Flag flag)
  : d_begin(NULL)
  , d_end(NULL)
  , d_gpu_data(NULL)
  , d_alloc_flag(flag)
{
    REQUIRE(count > 0);

#ifdef USE_CUDA
    unsigned int cuda_flags = 0;
    switch (d_alloc_flag)
    {
        case alloc::DEFAULT:
            cuda_flags = cudaHostAllocDefault;
            break;
        case alloc::MAPPED:
            cuda_flags = cudaHostAllocMapped;
            break;
        case alloc::WRITE_COMBINED:
            cuda_flags = cudaHostAllocWriteCombined;
            break;
        case alloc::MAPPED_WRITE_COMBINED:
            cuda_flags = cudaHostAllocMapped | cudaHostAllocWriteCombined;
            break;
    }

    // Allocate memory
    try
    {
        CudaCall(cudaHostAlloc(&d_begin, count * sizeof(T), cuda_flags));
    }
    catch (const profugus::assertion& e)
    {
        log(profugus::WARNING)
            << "Error: Failed to allocate "
            << count * sizeof(T) << " bytes on host";
        throw e;
    }

#else // CUDA is disabled in this build
    d_begin = reinterpret_cast<T*>(std::calloc(count, sizeof(T)));
#endif // USE_CUDA

    d_end = d_begin + count;
    std::fill(d_begin, d_end, value);

    if (is_mapped())
    {
#ifdef USE_CUDA
        // Get mapped GPU pointer for later access
        // (If mapped memory is unsupported on this system, this way we get an
        // error during construction rather than the kernel call.)
        CudaCall(cudaHostGetDevicePointer(&d_gpu_data, d_begin, 0));
#else // CUDA is disabled in this build
        std::free(d_begin);
        d_begin = NULL;
        INSIST(false,
                "Mapped memory cannot be emulated without CUDA support.");
#endif
    }

    ENSURE(count == size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Free memory on destruction
 */
template<class T>
Host_Vector<T>::~Host_Vector()
{
    try
    {
#ifdef USE_CUDA
        CudaCall(cudaFreeHost(d_begin));
#else
        std::free(d_begin);
#endif
    }
    catch (const profugus::assertion& e)
    {
        std::cerr << "!!! Error freeing device data "
            << "at " << d_begin << ": " << e.what() << std::endl;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory
 */
template<class T>
void Host_Vector<T>::assign(const_View_Field_t hostvec)
{
    REQUIRE(hostvec.size() == size());

    std::memcpy(d_begin, &hostvec[0], size() * sizeof(T));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory from another Host vector
 */
template<class T>
void Host_Vector<T>::assign(const This& rhs)
{
    REQUIRE(rhs.size() == size());

    std::memcpy(d_begin, rhs.d_begin, size() * sizeof(T));
}

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // cuda_utils_Host_Vector_t_hh

//---------------------------------------------------------------------------//
//                   end of cuda_utils/Host_Vector.t.hh
//---------------------------------------------------------------------------//
