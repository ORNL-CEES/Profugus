//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Transporter.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Source_Transporter template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "Source_Transporter.t.cuh"
#include "Uniform_Source.t.cuh"
#include "Fission_Source.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry Mesh_Geom;

// Instantiate class on geometry types
template class Source_Transporter<Mesh_Geom>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.pt.cu
//---------------------------------------------------------------------------//
