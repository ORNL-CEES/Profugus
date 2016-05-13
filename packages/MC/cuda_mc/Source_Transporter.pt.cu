//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source_Transporter.pt.cc
 * \author Stuart Slattery
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Source_Transporter template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source_Transporter.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Source_Transporter<Mesh_Geometry>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.pt.cc
//---------------------------------------------------------------------------//
