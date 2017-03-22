//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_geometry/Cartesian_Mesh.cu
 * \author Steven Hamilton
 * \date   Mon Jul 21 16:55:00 2014
 * \brief  Cartesian_Mesh member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "Cartesian_Mesh.hh"
#include "utils/Container_Props.hh"
#include "utils/View_Field.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize cartesian mesh with cell edges.
 */
Cartesian_Mesh_DMM::Cartesian_Mesh_DMM(const Vec_Dbl& x_edges,
                                       const Vec_Dbl& y_edges,
                                       const Vec_Dbl& z_edges)
    : d_x_edges(x_edges)
    , d_y_edges(y_edges)
    , d_z_edges(z_edges)
{
    // Hard-coded to 3D
    d_dimension = 3;

    INSIST(profugus::is_sorted(x_edges.begin(), x_edges.end()),
           "Mesh along x axis is not monotonically increasing.");
    INSIST(profugus::is_sorted(y_edges.begin(), y_edges.end()),
           "Mesh along y axis is not monotonically increasing.");
    INSIST(profugus::is_sorted(z_edges.begin(), z_edges.end()),
           "Mesh along z axis is not monotonically increasing.");

    // Compute cell volumes on host
    int num_cells_x = x_edges.size()-1;
    int num_cells_y = y_edges.size()-1;
    int num_cells_z = z_edges.size()-1;
    int num_cells = num_cells_x * num_cells_y * num_cells_z;
    d_volumes.resize(num_cells,0.0);
    for( int cell_k = 0; cell_k < num_cells_z; ++cell_k )
    {
        double width_k = z_edges[cell_k+1] - z_edges[cell_k];
        REQUIRE( width_k > 0.0 );
        for( int cell_j = 0; cell_j < num_cells_y; ++cell_j )
        {
            double width_j = y_edges[cell_j+1] - y_edges[cell_j];
            REQUIRE( width_j > 0.0 );
            for( int cell_i = 0; cell_i < num_cells_x; ++cell_i )
            {
                double width_i = x_edges[cell_i+1] - x_edges[cell_i];
                REQUIRE( width_j > 0.0 );
                cell_type cell = cell_i + num_cells_x *
                    (cell_j + num_cells_y * cell_k);
                d_volumes[cell] = width_i * width_j * width_k;
            }
        }
    }

    // Compute bounding box
    using def::I; using def::J; using def::K;
    d_lower[I] = x_edges.front();
    d_lower[J] = y_edges.front();
    d_lower[K] = z_edges.front();
    d_upper[I] = x_edges.back();
    d_upper[J] = y_edges.back();
    d_upper[K] = z_edges.back();
}

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.cu
//---------------------------------------------------------------------------//
