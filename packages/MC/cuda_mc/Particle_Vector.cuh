//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Particle_Vector.cuh
 * \author Steven Hamilton
 * \date   Tue Mar 07 15:48:29 2017
 * \brief  Particle_Vector definition.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Particle_Vector_cuh
#define MC_cuda_mc_Particle_Vector_cuh

#define USE_AOS 1

#if USE_AOS

#include "Particle_Vector_AOS.cuh"

namespace cuda_mc
{
template <class Geom>
using Particle_Vector = Particle_Vector_AOS<Geom>;

template <class Geom>
using Particle_Vector_DMM = Particle_Vector_AOS_DMM<Geom>;

} // end namespace mc

#else

#endif


//---------------------------------------------------------------------------//
#endif // MC_cuda_mc_Particle_Vector_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Particle_Vector.cuh
//---------------------------------------------------------------------------//
