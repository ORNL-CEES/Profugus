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

#ifdef __CUDACC__
#include <curand_kernel.h>
#endif

#include "cuda_utils/Definitions.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//

#ifdef __CUDACC__
typedef curandState_t RNG_State_t;
#endif

//! Fission site structure for storing fission sites in k-code.
struct Fission_Site
{
    int                         m;
    cuda_profugus::Space_Vector r;
};


//---------------------------------------------------------------------------//

} // end namespace cuda_mc

#endif // cuda_mc_Definitions_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Definitions.hh
//---------------------------------------------------------------------------//
