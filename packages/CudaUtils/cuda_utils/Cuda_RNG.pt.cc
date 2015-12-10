//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Cuda_RNG.pt.cc
 * \author Stuart Slattery
 * \date   Thu Dec 10 09:42:28 2015
 * \brief  Cuda RNG host instantiations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "Pseudo_Cuda.hh"
#include "Cuda_RNG.hh"
#include "Launch_Args.t.hh"

namespace cuda
{

template class Cuda_RNG<cuda::arch::Host>;
template class Cuda_RNG_Vector_Fill<cuda::arch::Host>;

} // end namespace cuda

//---------------------------------------------------------------------------//
//              end of cuda_utils/Cuda RNG.pt.cc
//---------------------------------------------------------------------------//
