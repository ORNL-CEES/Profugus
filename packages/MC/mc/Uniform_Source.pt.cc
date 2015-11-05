//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Uniform_Source.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Uniform_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Uniform_Source<Core>;
//template class Uniform_Source<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.pt.cc
//---------------------------------------------------------------------------//
