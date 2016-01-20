//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Bank.pt.cu
 * \author Stuart Slattery
 * \brief  Bank template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Bank.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_profugus
{

template class Bank<Mesh_Geometry>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Bank.pt.cc
//---------------------------------------------------------------------------//
