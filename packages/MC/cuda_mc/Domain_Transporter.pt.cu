//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter.pt.cc
 * \author Stuart Slattery
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Domain_Transporter template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Domain_Transporter.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Domain_Transporter<Mesh_Geometry>;

} // end namespace cuda_profuguss

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.pt.cc
//---------------------------------------------------------------------------//
