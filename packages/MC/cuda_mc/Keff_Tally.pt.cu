//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Keff_Tally template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Keff_Tally.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{

template class Keff_Tally_DMM<cuda_profugus::Mesh_Geometry>;
template class Keff_Tally_DMM<cuda_profugus::Core>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.pt.cu
//---------------------------------------------------------------------------//
