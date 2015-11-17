// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Stream.cc
 * \author Seth R Johnson
 * \date   Tue Oct 01 14:37:01 2013
 * \brief  Stream member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Stream.hh"

#include <iostream>
#include "harness/DBC.hh"
#include "CudaDBC.hh"

using cuda::arch::Host;
using cuda::arch::Device;

namespace cuda
{
//---------------------------------------------------------------------------//
/*!
 * \brief Create the Stream
 */
template<class Arch_T>
Stream<Arch_T>::Stream()
    : d_use_count(NULL)
{
    create_impl();

    d_use_count = new long;
    *d_use_count = 1;

    ENSURE(d_use_count);
    ENSURE(use_count() > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy the Stream, incrementing the reference counter
 */
template<class Arch_T>
Stream<Arch_T>::Stream(const This& rhs)
    : d_stream(rhs.d_stream)
    , d_use_count(rhs.d_use_count)
{
    REQUIRE(d_use_count);

    // Increment the reference count
    ++(*d_use_count);

    ENSURE(use_count() > 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment using copy-and-swap idiom
 *
 * This ensures exception safety and reduces duplicate code. No deleting or
 * checking for whether anything's been assigned.
 */
template<class Arch_T>
Stream<Arch_T>& Stream<Arch_T>::operator=(const This& rhs)
{
    This temp(rhs);
    this->swap(temp);

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap references and stream representationss with another stream
 */
template<class Arch_T>
void Stream<Arch_T>::swap(This& rhs)
{
    std::swap(d_stream   , rhs.d_stream);
    std::swap(d_use_count, rhs.d_use_count);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destroy the stream
 */
template<class Arch_T>
Stream<Arch_T>::~Stream()
{
    REQUIRE(use_count() > 0);

    // Decrement the reference count
    --(*d_use_count);

    if (use_count() == 0)
    {
        destroy_impl();
        delete d_use_count;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait for the Stream to complete
 */
template<class Arch_T>
void Stream<Arch_T>::synchronize()
{
    synchronize_impl();
}

//---------------------------------------------------------------------------//
/*!
 * \brief See if all prior events in the stream have completed
 */
template<class Arch_T>
bool Stream<Arch_T>::is_complete() const
{
    return is_complete_impl();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of objects using this stream
 */
template<class Arch_T>
long Stream<Arch_T>::use_count() const
{
    ENSURE(d_use_count);
    return *d_use_count;
}

#ifdef USE_CUDA
//---------------------------------------------------------------------------//
// DEVICE CODE
//---------------------------------------------------------------------------//
//! Create stream
template<>
void Stream<Device>::create_impl()
{
    CudaCall(cudaStreamCreate(&d_stream));
}

//---------------------------------------------------------------------------//
//! Wait for stream to complete
template<>
void Stream<Device>::synchronize_impl()
{
    CudaCall(cudaStreamSynchronize(d_stream));
}

//---------------------------------------------------------------------------//
//! Whether all prior events have completed
template<>
bool Stream<Device>::is_complete_impl() const
{
    return cudaStreamQuery(d_stream) == cudaSuccess;
}

//---------------------------------------------------------------------------//
//! Destroy stream if the use count is zero
template<>
void Stream<Device>::destroy_impl()
{
    // Don't ever throw in a destructor!
    try
    {
        CudaCall(cudaStreamDestroy(d_stream));
    }
    catch (const profugus::assertion& e)
    {
        std::cerr << "!!! Error cleaning up Cuda_Stream handle "
            << "at " << d_stream << ": " << e.what() << std::endl;
    }
}

//---------------------------------------------------------------------------//
// INSTANTIATE
//---------------------------------------------------------------------------//
template class Stream<Device>;
#endif //USE_CUDA

//---------------------------------------------------------------------------//
// HOST CODE (mostly null-ops)
//---------------------------------------------------------------------------//
//! Create stream
template<>
void Stream<Host>::create_impl()
{
    /* * */
}

//---------------------------------------------------------------------------//
//! Wait for stream to complete
template<>
void Stream<Host>::synchronize_impl()
{
    /* * */
}

//---------------------------------------------------------------------------//
//! Whether all prior events have completed
template<>
bool Stream<Host>::is_complete_impl() const
{
    return true;
}

//---------------------------------------------------------------------------//
//! Destroy stream if the use count is zero
template<>
void Stream<Host>::destroy_impl()
{
    /* * */
}

//---------------------------------------------------------------------------//
// INSTANTIATE
//---------------------------------------------------------------------------//
template class Stream<Host>;

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of Stream.cc
//---------------------------------------------------------------------------//
