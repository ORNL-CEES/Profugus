//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Source.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fission_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fission_Source.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Fission_Source<Core>;
//template class Fission_Source<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Source.pt.cc
//---------------------------------------------------------------------------//
