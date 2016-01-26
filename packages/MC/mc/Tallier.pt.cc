//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tallier.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Tallier template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Tallier.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Tallier<Core>;
template class Tallier<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Tallier.pt.cc
//---------------------------------------------------------------------------//
