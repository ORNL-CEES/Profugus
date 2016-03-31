//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/KDE_Kernel.pt.cc
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  KDE_Kernel template instantiations
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KDE_Kernel.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class KDE_Kernel<Core>;
template class KDE_Kernel<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of KDE_Kernel.pt.cc
//---------------------------------------------------------------------------//
