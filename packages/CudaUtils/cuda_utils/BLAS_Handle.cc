//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/BLAS_Handle.cc
 * \author Seth R Johnson
 * \date   Wed Dec 11 20:35:41 2013
 * \brief  BLAS_Handle member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "BLAS_Handle.hh"

#include <iostream>
#include "BLAS_device_common.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// STATIC DATA
//---------------------------------------------------------------------------//

//! CUDA BLAS handle (global because we reuse it)
cublasHandle_t BLAS_Handle::d_global_handle = NULL;

//! Reference counter for BLAS instances
unsigned int   BLAS_Handle::d_handle_count = 0;

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor for CUDA BLAS handle
 *
 * If this is the first extant instance, it will call cublasCreate to
 * initialize the handle.
 */
BLAS_Handle::BLAS_Handle()
{
    if (d_handle_count == 0)
    {
        REQUIRE(d_global_handle == NULL);

        CudaBlasCall(cublasCreate(&d_global_handle));
    }
    ++d_handle_count;
    ENSURE(d_global_handle != NULL);
    ENSURE(d_handle_count > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for CUDA BLAS handle
 */
BLAS_Handle::BLAS_Handle(const This& rhs)
{
    REQUIRE(d_handle_count > 0);

    CHECK(d_global_handle != NULL);

    ++d_handle_count;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destroy handle if it's the last one
 */
BLAS_Handle::~BLAS_Handle()
{
    REQUIRE(d_handle_count > 0);
    --d_handle_count;

    if (d_handle_count == 0)
    {
        // Don't ever throw in a destructor!
        try
        {
            CudaBlasCall(cublasDestroy(d_global_handle));
            // Only clear the pointer if destruct was successful
            d_global_handle = NULL;
        }
        catch (const profugus::assertion& e)
        {
            std::cerr << "!!! Error cleaning up CUDA handle: "
                << e.what() << std::endl;
        }
    }
}

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of BLAS_Handle.cc
//---------------------------------------------------------------------------//
