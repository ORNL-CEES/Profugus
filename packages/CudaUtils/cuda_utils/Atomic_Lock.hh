//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Atomic_Lock.hh
 * \author Seth R Johnson
 * \date   Thu Aug 15 09:32:04 2013
 * \brief  Atomic_Lock class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Profugus_cuda_utils_Atomic_Lock_hh
#define Profugus_cuda_utils_Atomic_Lock_hh

#include "Atomic_Lock_Kernel.cuh"

namespace cuda
{

//===========================================================================//
/*!
 * \class Atomic_Lock
 * \brief Host-side data management for atomic locks
 *
 * If your kernel needs an atomic lock, create an Atomic_Lock object in your
 * host code (kernel data). Then, pass the result of the data() function to the
 * GPU: the result is a thin int pointer with methods for acquiring and
 * releasing the lock.
 */
//===========================================================================//
template<class Arch_Switch>
class Atomic_Lock
{
  public:
    //@{
    //! Arch-dependent typedefs
    typedef Arch_Switch                Arch_t;
    typedef Atomic_Lock_Kernel<Arch_t> Atomic_Lock_Device_t;
    //@}

  public:
    // Create and initialize memory if appropriate
    Atomic_Lock();

    // Clean up memory if appropriate
    ~Atomic_Lock();

  public:
    // >>> Device access

    //! Create a "Lock_Device" device struct that references our data.
    Atomic_Lock_Device_t data()
    {
        return Atomic_Lock_Device_t(d_data);
    }

  private:
    // GPU-allocated memory if this uses the Device switch.
    int* d_data;

  private:
    // Prevent copies
    Atomic_Lock(const Atomic_Lock&);
};

} // end namespace cuda

#endif // Profugus_cuda_utils_Atomic_Lock_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Atomic_Lock.hh
//---------------------------------------------------------------------------//
