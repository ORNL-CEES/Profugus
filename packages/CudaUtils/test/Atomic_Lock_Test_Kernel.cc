//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/Atomic_Lock_Test_Kernel.cc
 * \author Seth R Johnson
 * \date   Thu Aug 15 08:18:39 2013
 * \brief  Atomic_Lock_Test_Kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

// Include fake cuda runtime
#include "../cuda_utils/Pseudo_Cuda.hh"

// Generate the polyglot Lock kernel for host code
#include "Atomic_Lock_Test_Kernel.cu"

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Lock_Test_Kernel.cc
//---------------------------------------------------------------------------//
