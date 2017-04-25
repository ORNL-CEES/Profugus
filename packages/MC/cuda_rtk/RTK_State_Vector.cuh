//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_State_Vector.cuh
 * \author Steven Hamilton
 * \date   Tue Apr 11 14:46:09 2017
 * \brief  RTK_State_Vector class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_State_Vector_cuh
#define MC_cuda_rtk_RTK_State_Vector_cuh

#define USE_RTK_AOS 0

#if USE_RTK_AOS

#include "RTK_State_Vector_AOS.cuh"

namespace cuda_profugus
{

using RTK_State_Vector     = RTK_State_Vector_AOS;
using RTK_State_Vector_DMM = RTK_State_Vector_AOS_DMM;

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

#else

#include "RTK_State_Vector_SOA.cuh"

namespace cuda_profugus
{

using RTK_State_Vector     = RTK_State_Vector_SOA;
using RTK_State_Vector_DMM = RTK_State_Vector_SOA_DMM;

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

#endif // USE_RTK_AOS

//---------------------------------------------------------------------------//
#endif // MC_cuda_rtk_RTK_State_Vector_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_State_Vector.cuh
//---------------------------------------------------------------------------//
