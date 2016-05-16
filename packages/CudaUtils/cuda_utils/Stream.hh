//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Stream.hh
 * \author Seth R Johnson
 * \date   Tue Oct 01 14:37:01 2013
 * \brief  Stream class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Stream_hh
#define CudaUtils_cuda_utils_Stream_hh

#include "config.h"
#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
#include <cuda_runtime.h>
#endif // USE_CUDA && !PSEUDO_CUDA

#include "Definitions.hh"

namespace cuda
{
// Declare Event class
template <typename Arch_T> class Event;

//===========================================================================//
/*!
 * \class Stream
 * \brief Synchronization management for CUDA.
 *
 * Streams can be used to asynchronously launch kernels, copy data, perform
 * timing, etc. There will be overhead with creating one, and with destroying
 * it, so try to make it in one place and keep references to it.
 */
/*!
 * \example cuda_utils/test/tstStream.cc
 *
 * Test of Stream.
 */
//===========================================================================//

template <typename Arch_T>
class Stream
{
    typedef Stream<Arch_T> This;
  public:
    //! Architecture type
    typedef Arch_T Arch_t;

#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
    //! Stream type
    typedef cudaStream_t stream_t;
#else
    typedef int stream_t;
#endif

  private:
    // >>> DATA

    // Underlying CUDA handle
    stream_t d_stream;

    // Reference-counted integer
    long *d_use_count;

  public:
    // >>> MANAGEMENT

    // Create the Stream, with optional flags
    Stream();

    // Copy the stream, incrementing the reference counter
    Stream(const This& rhs);

    // Assign the stream
    This& operator=(const This& rhs);

    // Swap with another stream
    void swap(This& rhs);

    // Destroy the Stream
    ~Stream();

    // >>> FUNCTIONS

    // Wait for the Stream to complete
    void synchronize();

    // >>> ACCESSORS

    // See if we're complete
    bool is_complete() const;

    //! Accessor to get the stream
    stream_t handle() const { return d_stream; }

    // Number of objects using this stream
    long use_count() const;

  private:
    // >>> Implementation differences between device/host
    void create_impl();
    void synchronize_impl();
    bool is_complete_impl() const;
    void destroy_impl();
};

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two streams
 *
 * This will be useful for having a "page flipping"-type buffering system.
 */
template <typename Arch_T>
inline void swap(Stream<Arch_T>& left, Stream<Arch_T>& right)
{
    left.swap(right);
}

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // CudaUtils_cuda_utils_Stream_hh

//---------------------------------------------------------------------------//
//                 end of Stream.hh
//---------------------------------------------------------------------------//
