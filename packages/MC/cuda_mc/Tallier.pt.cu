//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Tallier template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Tallier.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

template class Tallier<cuda_profugus::Mesh_Geometry>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Tallier.pt.cu
//---------------------------------------------------------------------------//
