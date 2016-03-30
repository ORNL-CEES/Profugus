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
Mesh_Geometry::Mesh_Geometry(const Vec_Dbl &x_edges,
                             const Vec_Dbl &y_edges,
                             const Vec_Dbl &z_edges)
    : d_mesh(x_edges,y_edges,z_edges)
{
    d_lower_bounds.x = x_edges.front();
    d_upper_bounds.x = x_edges.back();

    d_lower_bounds.y = y_edges.front();
    d_upper_bounds.y = y_edges.back();

    d_lower_bounds.z = z_edges.front();
    d_upper_bounds.z = z_edges.back();
}

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/Mesh_Geometry.cu
//---------------------------------------------------------------------------//
