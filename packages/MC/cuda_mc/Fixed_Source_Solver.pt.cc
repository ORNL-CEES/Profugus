//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fixed_Source_Solver.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fixed_Source_Solver template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fixed_Source_Solver.t.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Fixed_Source_Solver<Mesh_Geometry>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.pt.cc
//---------------------------------------------------------------------------//
