//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cuda_utils/test/Exp_Kernel.cc
 * \author Seth R Johnson
 * \date   Fri Aug 16 10:20:43 2013
 * \brief  Exp_Kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

// Include fake cuda runtime
#include "../Pseudo_Cuda.hh"

// Generate the polyglot Lock kernel for host code
#include "Exp_Kernel.cu"

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Exp_Kernel.cc
//---------------------------------------------------------------------------//
