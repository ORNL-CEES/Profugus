//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Definitions.hh
 * \author Steven Hamilton
 * \date   Friday April 25 16:46:37 2014
 * \brief  Monte Carlo Definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Definitions_hh
#define cuda_mc_Definitions_hh

#include <curand_kernel.h>

namespace cuda_mc
{

//---------------------------------------------------------------------------//

typedef curandState_t RNG_State_t;


//---------------------------------------------------------------------------//

} // end namespace cuda_mc

#endif // cuda_mc_Definitions_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Definitions.hh
//---------------------------------------------------------------------------//
