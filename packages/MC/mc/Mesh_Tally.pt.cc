//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Mesh_Tally.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Mesh_Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Mesh_Tally.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Mesh_Tally<Core>;
template class Mesh_Tally<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Mesh_Tally.pt.cc
//---------------------------------------------------------------------------//
