//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Current_Tally.pt.cc
 * \author Steven Hamilton
 * \date   Thu Apr 28 20:19:42 2016
 * \brief  Current_Tally member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Current_Tally.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Current_Tally<Core>;
template class Current_Tally<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Current_Tally.pt.cc
//---------------------------------------------------------------------------//
