//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fission_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fission_Source.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

template class Fission_Source<cuda_profugus::Mesh_Geometry>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Fission_Source.pt.cu
//---------------------------------------------------------------------------//
