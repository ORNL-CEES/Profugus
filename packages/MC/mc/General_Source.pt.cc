//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/General_Source.pt.cc
 * \author Steven Hamilton
 * \date   Mon Apr 04 20:38:12 2016
 * \brief  General_Source member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "General_Source.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class General_Source<Core>;
template class General_Source<Mesh_Geometry>;


} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of General_Source.pt.cc
//---------------------------------------------------------------------------//
