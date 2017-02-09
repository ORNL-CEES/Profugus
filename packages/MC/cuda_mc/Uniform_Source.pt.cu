//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Uniform_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class Uniform_Source_DMM<cuda_profugus::Mesh_Geometry>;
template class Uniform_Source_DMM<cuda_profugus::Core>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.pt.cu
//---------------------------------------------------------------------------//
