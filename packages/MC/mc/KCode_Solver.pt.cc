//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/KCode_Solver.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  KCode_Solver template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KCode_Solver.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class KCode_Solver<Core>;
//template class KCode_Solver<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of KCode_Solver.pt.cc
//---------------------------------------------------------------------------//
