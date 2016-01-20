//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle.pt.cu
 * \author Stuart Slattery
 * \brief  Particle template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Particle.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Particle<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Particle.pt.cc
//---------------------------------------------------------------------------//
