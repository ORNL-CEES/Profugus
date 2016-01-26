//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Rebalance.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fission_Rebalance template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fission_Rebalance.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Fission_Rebalance<Core>;
template class Fission_Rebalance<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Rebalance.pt.cc
//---------------------------------------------------------------------------//
