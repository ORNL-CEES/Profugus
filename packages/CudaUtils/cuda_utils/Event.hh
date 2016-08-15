//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Event.hh
 * \author Seth R Johnson
 * \date   Tue Jul 09 11:21:31 2013
 * \brief  Event class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Event_hh
#define cuda_utils_Event_hh

#include <config.h>
#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
#include <cuda_runtime.h>
#endif // USE_CUDA && !PSEUDO_CUDA

#include "Definitions.hh"

namespace cuda_utils
{
// Declare Stream class
template <typename Arch_T> class Stream;

//===========================================================================//
/*!
 * \class Event
 * \brief Wrap the CUDA event class, taking charge of create/destroy.
 *
 * This wraps the error checking/creation/destruction of the cudaEvent_t
 * struct. We throw a profugus::assertion if any of the operations fail.
 */
//===========================================================================//
template <typename Arch_T>
class Event
{
  public:
    //! Architecture type
    typedef Arch_T Arch_t;

  public:
    // Interface
    void record();
    void synchronize();
    float elapsed_time_since(const Event& start);

  private:
    // USE ONLY SPECIALIZATIONS: don't instantiate the generic class
    // This will fail if you don't have CUDA installed and try to use the
    // Device template
    Event();
};

//===========================================================================//
// HOST SPECIALIZATION (just use timer functionality)
//===========================================================================//
template<>
class Event<cuda_utils::arch::Host>
{
    typedef Event<cuda_utils::arch::Host> This;
  public:
    //! Architecture type
    typedef cuda_utils::arch::Host Arch_t;
    typedef Stream<Arch_t>   Stream_t;

  private:
    // Wall clock time at which record() is called
    double d_event_time;

  public:
    // Create the event, with optional flags
    Event()
      : d_event_time(-1.)
    {
        /* * */
    }

    // Record an event
    void record();

    // Record an event
    void record(Stream_t&) { record(); }

    // Wait for the event to complete (null-op for Host)
    void synchronize() { /* * */ }

    // See if we're complete
    bool is_complete() const { return true; }

    // Calculate elapsed time relative to another event in milliseconds
    float elapsed_time_since(const This& start) const;
};

#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
//===========================================================================//
// DEVICE SPECIALIZATION (actually use CUDA)
//===========================================================================//
template<>
class Event<cuda_utils::arch::Device>
{
    typedef Event<cuda_utils::arch::Device> This;
  public:
    //! Architecture type
    typedef cuda_utils::arch::Device Arch_t;
    typedef Stream<Arch_t>     Stream_t;

  private:
    // Underlying CUDA handle
    mutable cudaEvent_t d_handle;

  public:
    // Create the event, with optional flags
    explicit Event(unsigned int flags = cudaEventDefault);

    // Destroy the event
    ~Event();

    // Record an event (default stream)
    void record();

    // Record an event
    void record(Stream_t& stream);

    // Wait for the event to complete
    void synchronize();

    // See if we're complete
    bool is_complete() const;

    // Calculate elapsed time relative to another event in milliseconds
    float elapsed_time_since(const This& start) const;
};
#endif // USE_CUDA && !PSEUDO_CUDA

} // end namespace cuda_utils

#endif // cuda_utils_Event_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Event.hh
//---------------------------------------------------------------------------//
