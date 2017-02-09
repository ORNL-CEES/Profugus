//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Manager_Cuda.pt.cu
 * \author Steven Hamilton
 * \date   Wed Nov 25 11:42:29 2015
 * \brief  Manager_Cuda template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Manager_Cuda.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class Manager_Cuda<cuda_profugus::Mesh_Geometry_DMM>;
template class Manager_Cuda<cuda_profugus::Core_DMM>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Manager_cuda.pt.cu
//---------------------------------------------------------------------------//
