//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/BLAS_Handle.hh
 * \author Seth R Johnson
 * \date   Wed Dec 11 20:35:41 2013
 * \brief  BLAS_Handle class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * CUDA must be enabled to include this header.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_BLAS_Handle_hh
#define cuda_utils_BLAS_Handle_hh

#include <cublas_v2.h>

#include "harness/DBC.hh"

namespace cuda_utils
{

//===========================================================================//
/*!
 * \class BLAS_Handle
 * \brief Reference counter for a global CUDA handle.
 *
 * This should be used for internal purposes only. We don't include this in the
 * BLAS class because it needs to be global, not templated.
 */
//===========================================================================//
class BLAS_Handle
{
    typedef BLAS_Handle This;
  public:
    // Create the handle if the first instance
    BLAS_Handle();
    BLAS_Handle(const This& rhs);

    // Assignment is a null-op
    This& operator=(const This& rhs) { return *this; }

    // Decrement the handle count, destroy if necessary
    ~BLAS_Handle();

  public:
    // Access the handle, which under the hood is actually a pointer
    cublasHandle_t handle() const
    {
        CHECK(d_handle_count > 0);
        return d_global_handle;
    }

  public:
    //! Number of extant handles (primarily for testing)
    static unsigned int handle_count() { return d_handle_count; }

  private:
    // >>> PRIVATE STATIC DATA

    static cublasHandle_t d_global_handle;
    static unsigned int d_handle_count;
};

} // end namespace cuda_utils

#endif // cuda_utils_BLAS_Handle_hh

//---------------------------------------------------------------------------//
//                 end of BLAS_Handle.hh
//---------------------------------------------------------------------------//
