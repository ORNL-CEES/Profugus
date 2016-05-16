//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/mesh/Mesh.cc
 * \author Thomas M. Evans
 * \date   Wednesday February 12 0:20:19 2014
 * \brief  Mesh member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include "harness/Soft_Equivalence.hh"
#include "Mesh.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR AND DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for 3D or 2D mesh.
 *
 * The edges input should have size \c Nc_dim+1 where \c Nc_dim is the number
 * of cells in the \f$(x,y,z)\f$ dimension.
 *
 * An empty z_edges entry will construct a 2-D mesh.
 *
 * \param x_edges vector of cell edges in x in cm
 * \param y_edges vector of cell edges in y in cm
 * \param z_edges vector of cell edges in z in cm
 *
 * \param I_block i-index of the parallel partition
 * \param J_block j-index of the parallel partition
 * \param num_K_blocks number of pipelining blocks in k
 */
Mesh::Mesh(const def::Vec_Dbl &x_edges,
           const def::Vec_Dbl &y_edges,
           const def::Vec_Dbl &z_edges,
           size_type           I_block,
           size_type           J_block,
           size_type           num_K_blocks)
    : d_dimension(z_edges.empty() ? 2 : 3)
{
    REQUIRE(x_edges.size() >= 2);
    REQUIRE(y_edges.size() >= 2);
    REQUIRE(z_edges.size() >= 2 || z_edges.empty());
    REQUIRE(z_edges.empty() || num_K_blocks > 0);

    using def::I; using def::J; using def::K;

    // Set edges
    d_edges[I] = x_edges;
    d_edges[J] = y_edges;
    d_edges[K] = z_edges;

    // Set block index on I/J axis; num blocks along Z
    d_blocks[I] = I_block;
    d_blocks[J] = J_block;
    d_blocks[K] = num_K_blocks;

    // build widths, edges, block counts, etc
    this->build_mesh();

    ENSURE(d_num_cells == (x_edges.size() - 1)
                        * (y_edges.size() - 1)
                        * (d_dimension == 3 ? z_edges.size() - 1 : 1));
    ENSURE(d_num_vertices == x_edges.size()
                           * y_edges.size()
                           * (d_dimension == 3 ? z_edges.size() : 1));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for 2D mesh.
 *
 * The edges input should have size \c Nc_dim+1 where \c Nc_dim is the number
 * of cells in the \f$(x,y,z)\f$ dimension.
 *
 * An empty z_edges entry will construct a 2-D mesh.
 *
 * \param x_edges vector of cell edges in x in cm
 * \param y_edges vector of cell edges in y in cm
 *
 * \param I_block i-index of the parallel partition
 * \param J_block j-index of the parallel partition
 */
Mesh::Mesh(const def::Vec_Dbl &x_edges,
           const def::Vec_Dbl &y_edges,
           size_type           I_block,
           size_type           J_block)
    : d_dimension(2)
{
    REQUIRE(x_edges.size() >= 2);
    REQUIRE(y_edges.size() >= 2);

    using def::I; using def::J;

    // Set edges
    d_edges[I] = x_edges;
    d_edges[J] = y_edges;

    // Set block index on I/J axis
    d_blocks[I] = I_block;
    d_blocks[J] = J_block;

    // build widths, edges, block counts, etc
    this->build_mesh();

    ENSURE(d_num_cells    == (x_edges.size() - 1) * (y_edges.size() - 1));
    ENSURE(d_num_vertices == x_edges.size() * y_edges.size());
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Find the cell that contains a point.
 *
 * On return the \c ijk parameter is filled with the local \f$(i,j,k)\f$
 * coordinates of the cell if the point is in the mesh.
 *
 * For points on cell edges, a \c -n is returned where \c n is the index into
 * the cell edges array of the edge.  For example:
 * \verbatim
   |-------|-------|-------|
   |       |       |       |
   +   0   *   1   |   2   o
   |       |       |       |
   |-------|-------|-------|
   0       1       2       3
 * \endverbatim
 *
 * The \c "+" entry for \e i would return 0, for \c "*" the return is -1, and
 * for \c "o" it is -3.
 *
 * If this is a 2-D mesh, the edges along the K direction are -0.5 and +0.5 (to
 * get the 'per cm' normalization right), so the input k coordinate should be
 * zero.
 *
 * \param r Point to test
 * \param ijk Dim_Vector that is filled with the logical cell indices where
 * the point resides if it is in the mesh
 *
 * \return true if point is mesh; false otherwise
 */
bool Mesh::find_cell(const Space_Vector &r,
                     Dim_Vector         &ijk) const
{
    using def::I; using def::J; using def::K;
    REQUIRE(dimension() == 2 ? profugus::soft_equiv(r[K], 0.) : true);

    // first check to see if the point is in the mesh before doing a search
    // for the cell index
    if (r[I] < d_edges[I].front() || r[I] > d_edges[I].back())
        return false;
    if (r[J] < d_edges[J].front() || r[J] > d_edges[J].back())
        return false;
    if (r[K] < d_edges[K].front() || r[K] > d_edges[K].back())
        return false;

    // we are in the mesh so find the cell indices (unrolled loop)
    def::Vec_Dbl::const_iterator itr;

    // x
    itr = std::lower_bound(d_edges[I].begin(), d_edges[I].end(), r[I]);
    if (*itr != r[I])
    {
        ijk[I] = itr - d_edges[I].begin() - 1;
        ENSURE(ijk[I] >= 0 && ijk[I] < num_cells_dim(I));
    }
    else
    {
        ijk[I] = d_edges[I].begin() - itr;
    }

    // y
    itr    = std::lower_bound(d_edges[J].begin(), d_edges[J].end(), r[J]);
    if (*itr != r[J])
    {
        ijk[J] = itr - d_edges[J].begin() - 1;
        ENSURE(ijk[J] >= 0 && ijk[J] < num_cells_dim(J));
    }
    else
    {
        ijk[J] = d_edges[J].begin() - itr;
    }

    // z
    itr    = std::lower_bound(d_edges[K].begin(), d_edges[K].end(), r[K]);
    if (*itr != r[K])
    {
        ijk[K] = itr - d_edges[K].begin() - 1;
        ENSURE(ijk[K] >= 0 && ijk[K] < num_cells_dim(K));
    }
    else
    {
        ijk[K] = d_edges[K].begin() - itr;
    }

    return true;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Builds number of cells, cell widths, inverse widths
 */
void Mesh::build_mesh()
{
    using def::K;

    REQUIRE(d_dimension == 2 || d_dimension == 3);
    REQUIRE(d_dimension == 2 || d_blocks[K] > 0);
#ifdef REQUIRE_ON
    for (unsigned int d = 0; d < d_dimension; ++d)
        REQUIRE(!d_edges[d].empty());
#endif


    if (d_dimension == 2)
    {
        // Set Z axis to a unit length centered about zero
        d_edges[K].resize(2);
        d_edges[K][0] = -0.5;
        d_edges[K][1] =  0.5;

        // Set number of K blocks to one
        d_blocks[K] = 1;
    }

    // >>> Calculate number of cell in each dimension
    for (unsigned int d = 0; d < 3; ++d)
    {
        d_N[d]  = d_edges[d].size() - 1;
        d_Nb[d] = d_N[d];
    }

    // Z block count is divided by number of k blocks (1 if 2D)
    d_Nb[K] /= d_blocks[K];

    // >>> Calculate number of cells
    d_num_cells    = 1;
    d_num_vertices = 1;
    for (unsigned int d = 0; d < 3; ++d)
    {
        d_num_cells    *= d_N[d];
        d_num_vertices *= d_edges[d].size();
    }

    // In the 2D case, our two "edges" along the K axis really only count as
    // one vertex in the plane
    if (d_dimension == 2)
        d_num_vertices /= 2;

    // >>> Build cell widths and inverse widths
    for (unsigned int d = 0; d < 3; ++d)
    {
        d_width[d].resize(d_N[d]);
        d_inv_width[d].resize(d_N[d]);

        for (int l = 0; l < d_N[d]; ++l)
        {
            // cell width
            d_width[d][l] = d_edges[d][l+1] - d_edges[d][l];

            CHECK(d_width[d][l] > 0.0);

            // inverse width
            d_inv_width[d][l] = 1.0 / d_width[d][l];
        }
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Mesh.cc
//---------------------------------------------------------------------------//
