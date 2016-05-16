//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Definitions.hh
 * \author Seth R Johnson
 * \date   Wed Aug 14 11:52:04 2013
 * \brief  Common definitions for the CUDA utilities.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Definitions_hh
#define CudaUtils_cuda_utils_Definitions_hh

namespace cuda
{
//===========================================================================//

namespace arch
{
//---------------------------------------------------------------------------//
/*!
 * \struct Host
 * \brief  Switch to specify compilation on host (CPU) architecture
 */
struct Host
{
};

/*!
 * \struct Device
 * \brief  Switch to specify compilation on device (GPU) architecture
 */
struct Device
{
};
//---------------------------------------------------------------------------//
} // end namespace arch

//===========================================================================//
} // end namespace cuda

#endif // CudaUtils_cuda_utils_Definitions_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Definitions.hh
//---------------------------------------------------------------------------//
