//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.pt.cu
 * \author Stuart Slattery
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Tallier template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Tallier.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Tallier<Mesh_Geometry>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Tallier.pt.cc
//---------------------------------------------------------------------------//
