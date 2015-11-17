// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Atomic_Lock.cc
 * \author Seth R Johnson
 * \date   Thu Aug 15 09:32:04 2013
 * \brief  Atomic_Lock member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Atomic_Lock.hh"

#include <iostream>
#include <cstddef>
#include "config.h"
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#include "harness/DBC.hh"
#include "CudaDBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// DEFINITIONS FOR HOST LOCK SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Create, but don't do anything
 *
 * Locking on host code is meaningless
 */
template<>
Atomic_Lock<arch::Host>::Atomic_Lock()
  : d_data(NULL)
{
    /* * */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clean up memory if appropriate
 */
template<>
Atomic_Lock<arch::Host>::~Atomic_Lock()
{
    /* * */
}
//---------------------------------------------------------------------------//

#ifdef USE_CUDA
//---------------------------------------------------------------------------//
// DEFINITIONS FOR DEVICE LOCK SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize with one-element global memory set to 0
 */
template<>
Atomic_Lock<arch::Device>::Atomic_Lock()
  : d_data(NULL)
{
    CudaCall(cudaMalloc(&d_data, 1 * sizeof(int)));

    int zero_int = 0;
    CudaCall(cudaMemcpy(d_data, &zero_int, 1 * sizeof(int),
                cudaMemcpyHostToDevice));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clean up one-byte memory allocation when destroyed
 */
template<>
Atomic_Lock<arch::Device>::~Atomic_Lock()
{
    try
    {
        CudaCall(cudaFree(d_data));
    }
    catch (const profugus::assertion& e)
    {
        std::cerr << "!!! Error freeing device lock data "
            << "at " << d_data << ": " << e.what() << std::endl;
    }
}
//---------------------------------------------------------------------------//
#endif // USE_CUDA

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of Atomic_Lock.cc
//---------------------------------------------------------------------------//
