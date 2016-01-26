//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Problem_Builder.pt.cc
 * \author Steven Hamilton
 * \date   Wed Nov 25 11:50:17 2015
 * \brief  Problem_Builder template instantiations.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Problem_Builder.t.hh"

#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace mc
{

template class Problem_Builder<profugus::Core>;
template class Problem_Builder<profugus::Mesh_Geometry>;

} // end namespace mc

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.pt.cc
//---------------------------------------------------------------------------//
