//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc_driver/Geometry_Builder.pt.cc
 * \author Steven Hamilton
 * \date   Wed Nov 25 12:58:58 2015
 * \brief  Geometry_Builder member definitions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Geometry_Builder.t.hh"

#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace mc
{

template class Geometry_Builder<profugus::Core>;
template class Geometry_Builder<profugus::Mesh_Geometry>;

} // end namespace mc

//---------------------------------------------------------------------------//
//                 end of Geometry_Builder.pt.cc
//---------------------------------------------------------------------------//
