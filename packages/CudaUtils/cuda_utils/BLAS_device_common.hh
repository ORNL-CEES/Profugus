//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/BLAS_device_common.hh
 * \author Seth R Johnson
 * \date   Wed Dec 11 20:38:47 2013
 * \brief  BLAS_device_common kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_BLAS_device_common_hh
#define cuda_utils_BLAS_device_common_hh

#include <sstream>
#include "harness/DBC.hh"
#include "CudaDBC.hh"

namespace cuda_utils
{
//---------------------------------------------------------------------------//

//! Checked call for CUSA BLAS functions
#define CudaBlasCall(COND) \
    do \
    { \
        cublasStatus_t result = (COND); \
        if (result != CUBLAS_STATUS_SUCCESS) \
        { \
            std::ostringstream msg; \
            msg << "CUDA BLAS error " << result; \
            ::cuda_utils::toss_cuda_cookies(#COND, msg.str().c_str(), \
                    __FILE__, __LINE__); \
        } \
    } while (0);


//---------------------------------------------------------------------------//
} // end namespace cuda_utils
#endif // cuda_utils_BLAS_device_common_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/BLAS_device_common.hh
//---------------------------------------------------------------------------//
