//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/KCode_Solver.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  KCode_Solver template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KCode_Solver.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class KCode_Solver<Mesh_Geometry>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of KCode_Solver.pt.cc
//---------------------------------------------------------------------------//
