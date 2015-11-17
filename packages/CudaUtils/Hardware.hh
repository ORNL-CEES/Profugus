//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Hardware.hh
 * \author Seth R Johnson
 * \date   Tue Jul 09 15:43:18 2013
 * \brief  Device class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Hardware_hh
#define cuda_utils_Hardware_hh

#include <string>
#include "harness/DBC.hh"

#include "Definitions.hh"

namespace cuda
{

//===========================================================================//
/*!
 * \class Hardware
 * \brief C++ static wrapper class for abstracting GPU/CPU hardware layer
 *
 * This not only handles initialization but also stores useful memory/processor
 * capabilities of the device so that the program can intelligently choose how
 * to partition itself at runtime.
 *
 * If CUDA is disabled in this build, the acquire method of the Device
 * specialization will always raise an error.
 *
 * The default class provides sensible values; a specialization on Device
 * allows for CUDA interaction.
 *
 * \tparam Arch_T Architecture (Device or Host)
 */
//===========================================================================//
template<typename Arch_T>
class Hardware;

//! Specialization on Host
template<>
class Hardware<cuda::arch::Host>
{
  public:
    // >>> ACQUISITION

    //! Default/Host device is always valid
    static bool valid_device_exists() { return true; }

    //! Set up the device
    static void acquire() { /* * */ }

    //! Have we got a device
    static bool have_acquired() { return true; }

    //! Device name
    static std::string name() { return "CPU"; }

    // >>> HARDWARE

    //! Number of streaming multiprocessors
    static unsigned int num_multiprocessors() { return 1; }

    //! Number of CUDA cores in each multiprocessor
    static unsigned int num_cores_per_mp() { return 1; }

    //! Limit on number of threads per multiprocessor
    static unsigned int max_threads_per_mp() { return 1; }

    //! Number of threads per warp
    static int num_threads_per_warp() { return 1; }

    //! Amount of const memory (huge)
    static unsigned int const_memory() { return 1 << 30; }

    //! Amount of shared memory (huge)
    static unsigned int shared_memory_per_mp() { return 1 << 30; }

    //! Total number of cores
    static unsigned int num_total_cores() { return 1; }

    // >>> SYNCHRONIZATION

    //! ENSURE asynchronous operations are complete
    static void synchronize() { /* * */ }

  private:
    // Prevent instantiation
    Hardware();
};

//===========================================================================//
//! Specialization on Device
template<>
class Hardware<cuda::arch::Device>
{
  public:
    //! How to choose the device. We may add one that does e.g. node-based
    enum Acquire_Method
    {
        HIGHEST_FLOPS = 0,
        HIGHEST_COMPUTE,
        END_ACQUIRE_METHOD
    };

  private:
    // >>> DATA

    //! Device ID; -1 if unset
    static int d_device_id;

    //! Device name
    static std::string d_name;

    //! Number of streaming multiprocessors
    static unsigned int d_multiProcessorCount;

    //! Amount of constant memory available
    static unsigned int d_totalConstMem;

    //! Amount of shared memory available in each SM (divisible among blocks)
    static unsigned int d_sharedMemPerBlock;

    //! Total number of cores per multiprocessor (max physical thread capacity)
    static unsigned int d_cores_per_mp;

    //! Maximum number of threads per multiprocessor (warp scheduling)
    static unsigned int d_maxThreadsPerMultiProcessor;

  public:
    // >>> ACQUISITION

    // Is there a GPU in this system and are we compiling with CUDA?
    static bool valid_device_exists();

    // Get the best device
    static void acquire(Acquire_Method method = HIGHEST_FLOPS);

    //! Have we got a device yet?
    static bool have_acquired() { return d_device_id != -1; }

    //! Device name
    static std::string name() { return d_name; }

    // >>> HARDWARE

    //! Number of streaming multiprocessors
    static unsigned int num_multiprocessors()
    {
        REQUIRE(have_acquired());
        ENSURE(d_multiProcessorCount > 0);
        return d_multiProcessorCount;
    }

    //! Number of CUDA cores in each multiprocessor
    static unsigned int num_cores_per_mp()
    {
        REQUIRE(have_acquired());
        ENSURE(d_cores_per_mp > 0);
        return d_cores_per_mp;
    }

    //! Number of threads per warp
    static int num_threads_per_warp() { return 32; }

    //! Limit on number of threads per multiprocessor
    static unsigned int max_threads_per_mp()
    {
        REQUIRE(have_acquired());
        ENSURE(d_maxThreadsPerMultiProcessor > 0);
        return d_maxThreadsPerMultiProcessor;
    }

    //! Total constant memory on the device, in bytes
    static unsigned int const_memory()
    {
        REQUIRE(have_acquired());
        ENSURE(d_totalConstMem > 0);
        return d_totalConstMem;
    }

    //! Total shared memory per multiprocessor, in bytes
    static unsigned int shared_memory_per_mp()
    {
        REQUIRE(have_acquired());
        ENSURE(d_sharedMemPerBlock > 0);
        return d_sharedMemPerBlock;
    }

    //! Total number of cores
    static unsigned int num_total_cores()
    {
        REQUIRE(have_acquired());
        unsigned int result = d_multiProcessorCount * d_cores_per_mp;
        return result;
    }

    // >>> SYNCHRONIZATION

    // ENSURE asynchronous operations are complete
    static void synchronize();

  private:
    // Prevent instantiation
    Hardware();

  private:
    static unsigned int num_cores_per_mp(int major, int minor);
};

} // end namespace cuda

#endif // cuda_utils_Hardware_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Hardware.hh
//---------------------------------------------------------------------------//
