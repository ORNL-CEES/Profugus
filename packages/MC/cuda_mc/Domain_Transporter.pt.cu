//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Domain_Transporter template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Domain_Transporter.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

template class Domain_Transporter<cuda_profugus::Mesh_Geometry>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.pt.cu
//---------------------------------------------------------------------------//
