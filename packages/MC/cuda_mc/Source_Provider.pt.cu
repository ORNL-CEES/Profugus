//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Source_Provider.pt.cu
 * \author Steven Hamilton
 * \date   Wed Apr 13 08:38:06 2016
 * \brief  Source_Provider template instantiations.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source_Provider.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{
//---------------------------------------------------------------------------//

template class Source_Provider<cuda_profugus::Mesh_Geometry>;
template class Source_Provider<cuda_profugus::Core>;

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Source_Provider.pt.cu
//---------------------------------------------------------------------------//
