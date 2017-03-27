//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.pt.cu
 * \author Stuart Slattery
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Fission_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fission_Source.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_profugus
{

template class Fission_Source<Mesh_Geometry>;
template class Fission_Source<Core>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Source.pt.cu
//---------------------------------------------------------------------------//
