//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/VR_Roulette.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  VR_Roulette template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "VR_Roulette.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class VR_Roulette<Core>;
template class VR_Roulette<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.pt.cc
//---------------------------------------------------------------------------//
