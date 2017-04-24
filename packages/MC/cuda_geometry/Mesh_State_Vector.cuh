//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/Mesh_State_Vector.cuh
 * \author Steven Hamilton
 * \date   Tue Apr 11 14:46:09 2017
 * \brief  Mesh_State_Vector class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_Mesh_State_Vector_cuh
#define MC_cuda_geometry_Mesh_State_Vector_cuh

#define USE_AOS 0

#if USE_AOS

#include "Particle_Vector_AOS.cuh"

namespace cuda_profugus
{

using Mesh_State_Vector = Mesh_State_Vector_AOS;

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

#else

#include "Particle_Vector_SOA.cuh"

namespace cuda_profugus
{

using Mesh_State_Vector = Mesh_State_Vector_SOA;

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus


//---------------------------------------------------------------------------//
#endif // MC_cuda_geometry_Mesh_State_Vector_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/Mesh_State_Vector.cuh
//---------------------------------------------------------------------------//
