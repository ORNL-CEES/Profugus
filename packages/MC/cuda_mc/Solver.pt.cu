//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Solver.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Solver template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Solver.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_profugus
{

template class Solver<Mesh_Geometry>;
template class Solver<Core>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Solver.pt.cc
//---------------------------------------------------------------------------//
