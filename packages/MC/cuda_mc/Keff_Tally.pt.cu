//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally.pt.cu
 * \author Stuart Slattery
 * \brief  Keff_Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Keff_Tally.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Keff_Tally<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.pt.cc
//---------------------------------------------------------------------------//
