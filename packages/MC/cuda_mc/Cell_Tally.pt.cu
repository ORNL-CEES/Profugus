//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Cell_Tally.pt.cu
 * \author Stuart Slattery
 * \brief  Cell_Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Cell_Tally.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_profugus
{

template class Cell_Tally<Mesh_Geometry>;
template class Cell_Tally<Core>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Cell_Tally.pt.cc
//---------------------------------------------------------------------------//
