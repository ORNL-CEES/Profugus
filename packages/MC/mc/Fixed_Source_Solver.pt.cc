//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fixed_Source_Solver.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fixed_Source_Solver template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fixed_Source_Solver.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Fixed_Source_Solver<Core>;
template class Fixed_Source_Solver<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.pt.cc
//---------------------------------------------------------------------------//
