//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fixed_Source_Solver.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fixed_Source_Solver template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fixed_Source_Solver.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class Fixed_Source_Solver<cuda_profugus::Mesh_Geometry>;
template class Fixed_Source_Solver<cuda_profugus::Core>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.pt.cu
//---------------------------------------------------------------------------//
