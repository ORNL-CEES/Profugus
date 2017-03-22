//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/Mesh_Geometry.cu
 * \author Steven Hamilton
 * \date   Tue Dec 15 14:10:10 2015
 * \brief  Mesh_Geometry class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Mesh_Geometry.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
Mesh_Geometry_DMM::Mesh_Geometry_DMM(const Vec_Dbl &x_edges,
                                     const Vec_Dbl &y_edges,
                                     const Vec_Dbl &z_edges)
    : d_mesh(x_edges,y_edges,z_edges)
    , d_reflect(6,0)
{
}

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/Mesh_Geometry.cu
//---------------------------------------------------------------------------//
