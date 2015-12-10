//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Cuda_RNG.pt.cc
 * \author Stuart Slattery
 * \date   Thu Dec 10 09:42:28 2015
 * \brief  Cuda RNG device instantiations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Cuda_RNG.hh"

#include <curand_kernel.h>

namespace cuda
{

template class Cuda_RNG<cuda::arch::Device>;
template class Cuda_RNG_Vector_Fill<cuda::arch::Device>;

} // end namespace cuda

//---------------------------------------------------------------------------//
//              end of cuda_utils/Cuda RNG.cc
//---------------------------------------------------------------------------//
