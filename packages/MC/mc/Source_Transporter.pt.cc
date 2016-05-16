//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Source_Transporter.pt.cc
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Source_Transporter template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source_Transporter.t.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

namespace profugus
{

template class Source_Transporter<Core>;
template class Source_Transporter<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.pt.cc
//---------------------------------------------------------------------------//
