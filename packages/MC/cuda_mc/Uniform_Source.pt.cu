//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.pt.cu
 * \author Steven Hamilton
 * \date   Thu Nov 05 11:14:30 2015
 * \brief  Uniform_Source template instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source.t.cuh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry Mesh_Geom;

template class Uniform_Source<Mesh_Geom>;

// Instantiate get_particles free function
template thrust::device_vector<Particle<Mesh_Geom> > get_particles(
        cuda::Shared_Device_Ptr<Uniform_Source<Mesh_Geom>> &source,
        thrust::device_vector<cuda_mc::RNG_State_t>        &rngs);

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.pt.cu
//---------------------------------------------------------------------------//
