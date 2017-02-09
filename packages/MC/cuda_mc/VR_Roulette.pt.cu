//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/VR_Roulette.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  VR_Roulette template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "VR_Roulette.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class VR_Roulette<cuda_profugus::Mesh_Geometry>;
template class VR_Roulette<cuda_profugus::Core>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.pt.cu
//---------------------------------------------------------------------------//
