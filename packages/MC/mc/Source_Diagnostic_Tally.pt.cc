//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Diagnostic_Tally.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Source_Diagnostic_Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source_Diagnostic_Tally.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Source_Diagnostic_Tally<Core>;
//template class Source_Diagnostic_Tally<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source_Diagnostic_Tally.pt.cc
//---------------------------------------------------------------------------//
