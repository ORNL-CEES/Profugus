//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Cuda_RNG.cc
 * \author Stuart Slattery
 * \date   Thu Dec 10 09:42:28 2015
 * \brief  Cuda RNG host instantiations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "Pseudo_Cuda.hh"
#include "Cuda_RNG.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// Global CUDA RNG controller.
//---------------------------------------------------------------------------//
// Initial declaration of the global controller.
profugus::RNG_Control Cuda_Global_RNG_Control::d_rng_control(0);

// Initialize the global controller with a global seed.
void Cuda_Global_RNG_Control::initialize( const int global_seed )
{
    d_rng_control = 
	profugus::RNG_Control( global_seed, std::numeric_limits<int>::max() );
}

} // end namespace cuda

//---------------------------------------------------------------------------//
//              end of cuda_utils/Cuda RNG.cc
//---------------------------------------------------------------------------//
