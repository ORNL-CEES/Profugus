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
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class Tallier_DMM<cuda_profugus::Mesh_Geometry>;
template class Tallier_DMM<cuda_profugus::Core>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Tallier.pt.cu
//---------------------------------------------------------------------------//
