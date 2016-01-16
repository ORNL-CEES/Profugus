//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.pt.cu
 * \author Stuart Slattery
 * \brief  Uniform_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source.t.cuh"
#include "Box_Shape.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Uniform_Source<Mesh_Geometry,Box_Shape>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.pt.cc
//---------------------------------------------------------------------------//
