//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/Profiler_Kernel.cc
 * \author Seth R Johnson
 * \date   Sun Aug 18 07:13:45 2013
 * \brief  Profiler_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

// Include fake cuda runtime
#include "../cuda_utils/Pseudo_Cuda.hh"

// Generate the polyglot Lock kernel for host code
#include "Profiler_Kernel.cu"

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Profiler_Kernel.cc
//---------------------------------------------------------------------------//
