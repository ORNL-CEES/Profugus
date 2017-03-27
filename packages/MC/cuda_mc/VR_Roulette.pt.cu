//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/VR_Roulette.pt.cu
 * \author Stuart Slattery
 * \brief  VR_Roulette template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "VR_Roulette.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_profugus
{

template class VR_Roulette<Mesh_Geometry>;
template class VR_Roulette<Core>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.pt.cc
//---------------------------------------------------------------------------//
