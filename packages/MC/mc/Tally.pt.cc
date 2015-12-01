//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tally.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Tally.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Tally<Core>;
template class Tally<Mesh_Geometry>;

template class Source_Tally<Core>;
template class Source_Tally<Mesh_Geometry>;

template class Pathlength_Tally<Core>;
template class Pathlength_Tally<Mesh_Geometry>;

template class Compound_Tally<Core>;
template class Compound_Tally<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Tally.pt.cc
//---------------------------------------------------------------------------//
