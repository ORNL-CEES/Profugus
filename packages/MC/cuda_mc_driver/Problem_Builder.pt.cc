//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Problem_Builder.pt.cc
 * \author Steven Hamilton
 * \date   Wed Nov 25 11:50:17 2015
 * \brief  Problem_Builder template instantiations.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Problem_Builder.t.hh"

#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

template class Problem_Builder<cuda_profugus::Mesh_Geometry>;

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.pt.cc
//---------------------------------------------------------------------------//
