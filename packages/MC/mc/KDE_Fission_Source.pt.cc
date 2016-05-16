//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/KDE_Fission_Source.pt.cc
 * \author Gregory Davidson
 * \date   Mon Nov 23 15:47:23 2015
 * \brief  KDE_Fission_Source template instantiations.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KDE_Fission_Source.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class KDE_Fission_Source<Core>;
template class KDE_Fission_Source<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//              end of KDE_Fission_Source.pt.cc
//---------------------------------------------------------------------------//
