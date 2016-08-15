// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Event.cc
 * \author Seth R Johnson
 * \date   Tue Jul 09 11:21:31 2013
 * \brief  Event member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Event.hh"

#include <iostream>
#include "harness/DBC.hh"
#include "comm/Timer.hh"

#include "CudaDBC.hh"
#include "Hardware.hh"
#include "Stream.hh"

using cuda_utils::arch::Host;
using cuda_utils::arch::Device;

namespace cuda_utils
{
//===========================================================================//
// HOST SPECIALIZATION
//===========================================================================//
/*!
 * \brief Record an event
 */
void Event<Host>::record()
{
    // Note: These definitions and their headers should have been propagated
    // through Timer.hh
#if defined(HAVE_SYS_TIMES_H)
    tms tms_now;
    d_event_time = profugus::wall_clock_time(tms_now);
#elif defined(HAVE_WINDOWS_H)
    // Timer frequency and current time
    long long ll_freq, ll_now;
    QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&ll_freq));
    QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(&ll_now));
    CHECK(ll_freq > 0);
    d_event_time = static_cast<double>(ll_now) / ll_freq;
#else
    INSIST(false, "Timing is not available on this platform.");
#endif
    ENSURE(d_event_time > 0.);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate elapsed time relative to another event (ms)
 *
 * This should be the "stop" event.
 *
 * \return Elapsed time in milliseconds.
 */
float Event<Host>::elapsed_time_since(const This& start) const
{
    return (d_event_time - start.d_event_time) * 1000;
}

#ifdef USE_CUDA
//===========================================================================//
// DEVICE SPECIALIZATION
//===========================================================================//
/*!
 * \brief Create the event, with optional flags
 *
 * \pre Device must already have been set up.
 */
Event<Device>::Event(unsigned int flags)
 : d_handle(NULL)
{
    typedef cuda_utils::Hardware<Arch_t> Hardware_t;
    REQUIRE(Hardware_t::have_acquired());

    CudaCall(cudaEventCreateWithFlags(&d_handle, flags));
    ENSURE(d_handle != NULL);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destroy the event
 */
Event<Device>::~Event()
{
    if (d_handle != NULL)
    {
        // Don't ever throw in a destructor!
        try
        {
            CudaCall(cudaEventDestroy(d_handle));
        }
        catch (const profugus::assertion& e)
        {
            std::cerr << "!!! Error cleaning up Cuda_Event handle "
                << "at " << d_handle << ": " << e.what() << std::endl;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Record this event in the default stream
 */
void Event<Device>::record()
{
    CudaCall(cudaEventRecord(d_handle));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Record this event in a given stream
 */
void Event<Device>::record(Stream_t& stream)
{
    CudaCall(cudaEventRecord(d_handle, stream.handle()));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait for the event to complete
 */
void Event<Device>::synchronize()
{
    CudaCall(cudaEventSynchronize(d_handle));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the event is complete
 */
bool Event<Device>::is_complete() const
{
    return cudaEventQuery(d_handle) == cudaSuccess;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate elapsed time relative to another event (ms)
 *
 * This should be the "stop" event.
 *
 * \return Elapsed time in milliseconds.
 */
float Event<Device>::elapsed_time_since(const This& start) const
{
    float time;
    CudaCall(cudaEventElapsedTime(&time, start.d_handle, this->d_handle));
    return time;
}

#endif //USE_CUDA
//---------------------------------------------------------------------------//
} // end namespace cuda_utils

//---------------------------------------------------------------------------//
//                 end of Event.cc
//---------------------------------------------------------------------------//
