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
#include "cuda_utils/Host_Vector.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize cartesian mesh with cell edges.
 */
Cartesian_Mesh::Cartesian_Mesh(const Vec_Dbl& x_edges,
                               const Vec_Dbl& y_edges,
                               const Vec_Dbl& z_edges)
    : d_x_edges_vec(profugus::make_view(x_edges))
    , d_y_edges_vec(profugus::make_view(y_edges))
    , d_z_edges_vec(profugus::make_view(z_edges))
{
    // Hard-coded to 3D
    d_dimension = 3;

    // Get cell counts for each dimensions
    d_cells_x = x_edges.size() - 1;
    d_cells_y = y_edges.size() - 1;
    d_cells_z = z_edges.size() - 1;
    REQUIRE( d_cells_x > 0 );
    REQUIRE( d_cells_y > 0 );
    REQUIRE( d_cells_z > 0 );

    // Get raw pointers from device vectors
    dd_x_edges = d_x_edges_vec.data();
    dd_y_edges = d_y_edges_vec.data();
    dd_z_edges = d_z_edges_vec.data();

    INSIST(profugus::is_sorted(x_edges.begin(), x_edges.end()),
           "Mesh along x axis is not monotonically increasing.");
    INSIST(profugus::is_sorted(y_edges.begin(), y_edges.end()),
           "Mesh along y axis is not monotonically increasing.");
    INSIST(profugus::is_sorted(z_edges.begin(), z_edges.end()),
           "Mesh along z axis is not monotonically increasing.");

    d_num_cells = d_cells_x * d_cells_y * d_cells_z;

    // Compute cell volumes on host
    cuda::Host_Vector<double> host_volumes(d_num_cells);
    for( int cell_k = 0; cell_k < d_cells_z; ++cell_k )
    {
        double width_k = z_edges[cell_k+1] - z_edges[cell_k];
        REQUIRE( width_k > 0.0 );
        for( int cell_j = 0; cell_j < d_cells_y; ++cell_j )
        {
            double width_j = y_edges[cell_j+1] - y_edges[cell_j];
            REQUIRE( width_j > 0.0 );
            for( int cell_i = 0; cell_i < d_cells_x; ++cell_i )
            {
                double width_i = x_edges[cell_i+1] - x_edges[cell_i];
                REQUIRE( width_j > 0.0 );
                cell_type cell;
                this->index(cell_i,cell_j,cell_k,cell);
                host_volumes[cell] = width_i * width_j * width_k;
            }
        }
    }
    d_volumes_vec = std::make_shared<Device_Vector<double> >(d_num_cells);
    d_volumes_vec->assign(host_volumes);
    dd_volumes = d_volumes_vec->data();

    ENSURE(d_num_cells >= 1);
}

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.cu
//---------------------------------------------------------------------------//
