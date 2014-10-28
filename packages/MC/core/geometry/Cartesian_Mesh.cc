//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/geometry/Cartesian_Mesh.cc
 * \author Thomas M. Evans
 * \date   Mon Jul 21 16:55:00 2014
 * \brief  Cartesian_Mesh member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Cartesian_Mesh.hh"

#include "utils/Container_Functions.hh"

//---------------------------------------------------------------------------//
// ANONYMOUS NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//

namespace
{

//---------------------------------------------------------------------------//
/*!
 * \brief Small private function to return index with sign.
 *
 * If val is less than or equal to the first element, the result is -1; if
 * greater than the last element, the result is N.
 */
inline int edge_signed_index(const std::vector<double>& edges,
                             double                     val)
{
    std::vector<double>::const_iterator result =
        std::lower_bound(edges.begin(), edges.end(), val);
    CHECK(result != edges.end());

    if (*result != val)
    {
        // not on an edge: positive cell
        return result - edges.begin() - 1;
    }
    else
    {
        // not on an edge: negative cell
        return edges.begin() - result;
    }
}

} // end anonymous namespace

//---------------------------------------------------------------------------//
// PROFUGUS NAMESPACE
//---------------------------------------------------------------------------//

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize cartesian mesh with cell edges.
 *
 * To make this a 2-D mesh, leave the z edges empty.
 */
Cartesian_Mesh::Cartesian_Mesh(const Vec_Dbl& x_edges,
                               const Vec_Dbl& y_edges,
                               const Vec_Dbl& z_edges)
    : d_edges(x_edges, y_edges, z_edges)
{
    d_num_cells = 1;
    d_dimension = 0;

    for (int ax = 0; ax < 3; ++ax)
    {
        if (d_edges[ax].empty())
        {
            d_extents[ax] = 1;
            break;
        }
        // Set extents along this axis
        d_extents[ax] = d_edges[ax].size() - 1;
        // Set number of cells
        d_num_cells *= d_extents[ax];
        // Increase dimension
        ++d_dimension;

        VALIDATE(profugus::is_sorted(d_edges[ax].begin(), d_edges[ax].end()),
                 "Mesh along " << "xyz"[ax] << " axis is not monotonically "
                 "increasing.");
    }

    // For now only support 2-D (missing Z) or 3-D
    VALIDATE(d_dimension == 3 || (d_dimension == 2 && d_edges[def::Z].empty()),
             "Only xy and xyz meshes are currently supported.");

    ENSURE(d_dimension >= 2);
    ENSURE(d_num_cells >= 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a cardinal mesh index to an (i,j,k) index.
 *
 * \return cardinal mesh index.
 */
void Cartesian_Mesh::cardinal(size_type cell,
                              dim_type& i,
                              dim_type& j,
                              dim_type& k) const
{
    REQUIRE(cell < num_cells());

    using def::I; using def::J; using def::K;

    k = cell / (d_extents[I] * d_extents[J]);
    cell -= k * (d_extents[I] * d_extents[J]);

    j = cell / d_extents[I];
    cell -= j * d_extents[I];

    i = cell;

    ENSURE(0 <= i && i < d_extents[I]);
    ENSURE(0 <= j && j < d_extents[J]);
    ENSURE(0 <= k && k < d_extents[K]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Locate position's ijk coordinates with upper edges being "inside".
 *
 * If the spatial coordinate's position is less than or equal to the first
 * element, the result is -1; if greater than the last element, the result is
 * N (number of cells).
 *
 * This allows for testing of being inside the mesh or on what outside portion
 * of the mesh it is but it does not give information about being on an edge or
 * not.
 */
void Cartesian_Mesh::find_upper(const Space_Vector& r, Dim_Vector& ijk) const
{
    using def::I; using def::J; using def::K;

    ijk[I] = find_upper(r[I], I);
    ijk[J] = find_upper(r[J], J);
    ijk[K] = find_upper(r[K], K);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Locate position's ijk coordinates.
 *
 * If the position is on an edge, return the negative cell index. See kba mesh
 * for a description.
 */
bool Cartesian_Mesh::find(const Space_Vector& r, Dim_Vector& ijk) const
{
    using def::I; using def::J; using def::K;

    // Test for being inside the edges
    if (    r[I] < d_edges[I].front() || r[I] >= d_edges[I].back()
         || r[J] < d_edges[J].front() || r[J] >= d_edges[J].back()
         || r[K] < d_edges[K].front() || r[K] >= d_edges[K].back())
        return false;

    ijk[I] = edge_signed_index(d_edges[I], r[I]);
    ijk[J] = edge_signed_index(d_edges[J], r[J]);
    ijk[K] = edge_signed_index(d_edges[K], r[K]);
    return true;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.cc
//---------------------------------------------------------------------------//
