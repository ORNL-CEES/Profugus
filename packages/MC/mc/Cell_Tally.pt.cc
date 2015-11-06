//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Cell_Tally.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Cell_Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Cell_Tally.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Cell_Tally<Core>;
//template class Cell_Tally<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Cell_Tally.pt.cc
//---------------------------------------------------------------------------//
