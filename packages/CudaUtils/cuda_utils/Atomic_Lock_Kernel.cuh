//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Atomic_Lock_Kernel.cuh
 * \author Seth R Johnson
 * \date   Thu Aug 15 07:57:35 2013
 * \brief  Atomic_Lock_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Atomic_Lock_Kernel_cuh
#define cuda_utils_Atomic_Lock_Kernel_cuh

#include "Definitions.hh"

namespace cuda
{
// Declare class so it can be our friend
template<typename Arch_Switch> class Atomic_Lock;

//===========================================================================//
// CLASS DEFINITIONS
//===========================================================================//
/*!
 * \class Atomic_Lock_Kernel
 * \brief Device-wide synchronization
 *
 * This is used in kernel code to guarantee that two different threadblocks
 * aren't running a section of the code at the same time.
 *
 * The lock object must be created on the host, and then passed
 * to the kernel.
 *
 * \note We recommend that the \c __syncthreads command is used so that all
 * warps in a block reach the locking commands at the same time. (This allows
 * the SM to do other work before hitting the spin lock.)
 *
 * Example:
    \code
    Atomic_Lock_Kernel<Arch_Switch> lock;

    __syncthreads();
    if (threadIdx.x == 0)
        lock.acquire();
    __syncthreads();

    // Do something in *only* the threadblock that has the lock

    __syncthreads();
    if (threadIdx.x == 0)
        lock.release();
    __syncthreads();

    \endcode
 */
//===========================================================================//
template<class Arch_Switch>
class Atomic_Lock_Kernel
{
  private:
    // Only use specializations!
    Atomic_Lock_Kernel();
};

//===========================================================================//
/*!
 * \brief Specialization for device code
 *
 * This thin structure uses a raw pointer reference to device memory *which
 * should be allocated by Atomic_Lock*.  We don't allocate it, we don't
 * validate it, we don't delete it.
 */
template<>
class Atomic_Lock_Kernel<arch::Device>
{
    typedef Atomic_Lock_Kernel<arch::Device> This;
  private:
    // >>> DATA

    //! Device memory: whether we're locked across all blocks in this kernel.
    int* d_lock;

  private:
    // >>> CONSTRUCTION

    friend class Atomic_Lock<arch::Device>;

    // Initialize on the host (only our friend can create)
    Atomic_Lock_Kernel(int* const lock) : d_lock(lock) { /* * */ }

  public:
    //! Implicit copy constructor during kernel call
    Atomic_Lock_Kernel(const This& rhs) : d_lock(rhs.d_lock) { /* * */ }

  public:
    // >>> LOCKING

#ifdef __NVCC__
    //! Wait (spin) until no other threadblock is using the Atomic_Lock_Kernel
    __device__ void acquire()
    {
        // Locked = 1, unlocked = 0
        while (atomicCAS(d_lock, 0, 1))
        {
            /* * */
        }
    }

    //! Allow other threadblocks to continue
    __device__ void release()
    {
        *d_lock = 0;
    }
#endif
};

//===========================================================================//
//! Specialization for CPU emulation code
template<>
class Atomic_Lock_Kernel<arch::Host>
{
    typedef Atomic_Lock_Kernel<arch::Host> This;
  private:
    // >>> CONSTRUCTION

    friend class Atomic_Lock<arch::Host>;

    // Initialize on the host (only our friend can create)
    Atomic_Lock_Kernel(int* const)  { /* * */ }

  public:
    //! Implicit copy constructor during kernel call
    Atomic_Lock_Kernel(const This& rhs) { /* * */ }

    // >>> LOCKING

    //! In single-thread CPU code, no locking is needed
    void acquire() { /* * */ }

    //! In single-thread CPU code, no locking is needed
    void release() { /* * */ }
};

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_Atomic_Lock_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Atomic_Lock_Kernel.cuh
//---------------------------------------------------------------------------//
