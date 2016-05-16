//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/geometry/Cartesian_Mesh.i.hh
 * \author Thomas M. Evans
 * \date   Mon Jul 21 16:55:00 2014
 * \brief  Member definitions of class Cartesian_Mesh.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_geometry_Cartesian_Mesh_i_hh
#define MC_geometry_Cartesian_Mesh_i_hh

#include <algorithm>

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Convert global cell index to dimension vector
 */
Cartesian_Mesh::Dim_Vector Cartesian_Mesh::cardinal(size_type cell) const
{
    REQUIRE(cell < num_cells());
    Dim_Vector result;
    cardinal(cell, result[0], result[1], result[2]);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the volume of a cell
 */
double Cartesian_Mesh::volume(int i, int j, int k) const
{
    using def::I; using def::J; using def::K;

    REQUIRE(i >= 0 && i < d_extents[I]);
    REQUIRE(j >= 0 && j < d_extents[J]);
    REQUIRE(k >= 0 && k < d_extents[K]);

    return (d_edges[I][i+1] - d_edges[I][i])
         * (d_edges[J][j+1] - d_edges[J][j])
         * (d_edges[K][k+1] - d_edges[K][k]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute volume from global cell index (conveneience function)
 */
double Cartesian_Mesh::volume(size_type global_cell_index) const
{
    dim_type i, j, k;
    cardinal(global_cell_index, i, j, k);
    return this->volume(i, j, k);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a global cell index from dimensional indices
 *
 * This automatically performs the correct testing of -1 or N edges that
 * \c find_upper returns.
 *
 * \return success of finding cell
 */
bool Cartesian_Mesh::index(dim_type i, dim_type j, dim_type k,
        size_type& cell) const
{
    using def::I; using def::J; using def::K;

    // Check whether the dimesions are inside the mesh
    if (    i < 0 || i >= d_extents[I]
         || j < 0 || j >= d_extents[J]
         || k < 0 || k >= d_extents[K])
        return false;

    // Calculate cell index
    cell = i + d_extents[I] * (j + k * d_extents[J]);
    ENSURE(cell < num_cells());
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Locate cell index along an axis with upper edges being "inside"
 *
 * If the spatial coordinate's position is less than or equal to the first
 * element, the result is -1; if greater than the last element, the result is
 * N (number of cells).
 */
Cartesian_Mesh::dim_type
Cartesian_Mesh::find_upper(double r, dim_type d) const
{
    REQUIRE(0 <= d && d < d_dimension);
    return std::lower_bound(d_edges[d].begin(), d_edges[d].end(), r)
        - d_edges[d].begin() - 1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a global cell index from dimensional indices
 *
 * This function raises an error if the inputs are out of bounds.
 *
 * \return success of finding cell
 */
Cartesian_Mesh::size_type
Cartesian_Mesh::index(dim_type i, dim_type j, dim_type k) const
{
    using def::I; using def::J; using def::K;

    // Check whether the dimesions are inside the mesh
    REQUIRE(i >= 0 && i < d_extents[I]);
    REQUIRE(j >= 0 && j < d_extents[J]);
    REQUIRE(k >= 0 && k < d_extents[K]);

    // Calculate cell index
    return i + d_extents[I] * (j + k * d_extents[J]);
}

} // end namespace profugus

#endif // MC_geometry_Cartesian_Mesh_i_hh

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.i.hh
//---------------------------------------------------------------------------//
