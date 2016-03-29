//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Axial_KDE_Kernel.pt.cc
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  Axial_KDE_Kernel template instantiations
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Axial_KDE_Kernel.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Axial_KDE_Kernel<Core>;
template class Axial_KDE_Kernel<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Axial_KDE_Kernel.pt.cc
//---------------------------------------------------------------------------//
