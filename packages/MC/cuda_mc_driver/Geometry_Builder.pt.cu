//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Geometry_Builder.pt.cu
 * \author Steven Hamilton
 * \date   Wed Nov 25 12:58:58 2015
 * \brief  Geometry_Builder member definitions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Geometry_Builder.t.cuh"

#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class Geometry_Builder<cuda_profugus::Mesh_Geometry_DMM>;
template class Geometry_Builder<cuda_profugus::Core_DMM>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Geometry_Builder.pt.cu
//---------------------------------------------------------------------------//
