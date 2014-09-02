//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Mesh.i.hh
 * \author Thomas M. Evans
 * \date   Wednesday February 12 0:21:50 2014
 * \brief  Member definitions of class Mesh.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Mesh_i_hh
#define spn_Mesh_i_hh

#include "harness/Soft_Equivalence.hh"

namespace profugus
{


//---------------------------------------------------------------------------//
// PUBLIC INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Convert an (i,j) index to a cardinal mesh index assuming k=0.
 *
 * \return cardinal mesh index.
 */
Mesh::size_type Mesh::convert(size_type i,
                              size_type j) const
{
    REQUIRE(static_cast<int>(i) < num_cells_dim(def::I) );
    REQUIRE(static_cast<int>(j) < num_cells_dim(def::J) );

    return i + j * num_cells_dim(def::I);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert an (i,j,k) index to a cardinal mesh index.
 *
 * \return cardinal mesh index.
 */
Mesh::size_type Mesh::convert(size_type i,
                              size_type j,
                              size_type k) const
{
    REQUIRE(static_cast<int>(i) < num_cells_dim(def::I) );
    REQUIRE(static_cast<int>(j) < num_cells_dim(def::J) );

    // Allow 3 arguments for 2D mesh if k=0
    REQUIRE( d_dimension==3 ? static_cast<int>(k) < num_cells_dim(def::K) :
              k == 0 );

    return i + j * num_cells_dim(def::I)
        + k * num_cells_dim(def::I) * num_cells_dim(def::J);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a cardinal mesh index to an (i,j,k) index.
 *
 * \return cardinal mesh index.
 */
Mesh::Dim_Vector Mesh::cardinal(size_type cell) const
{
    REQUIRE(cell < num_cells() );

    Dim_Vector ijk;
    size_type index = cell;
    ijk[def::K] = index /
                  ( num_cells_dim(def::I)*num_cells_dim(def::J) );
    index      -= ijk[def::K] *
                  ( num_cells_dim(def::I)*num_cells_dim(def::J) );
    ijk[def::J] = index / num_cells_dim(def::I);
    index      -= ijk[def::J] * num_cells_dim(def::I);
    ijk[def::I] = index;
    return ijk;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the volume of a cell by cellid.
 *
 * \param cell index
 *
 * \return volume of the cell in cm^3
 */
double Mesh::volume(size_type cell) const
{
    REQUIRE(cell < num_cells());

    Dim_Vector ijk = cardinal(cell);

    return (width(ijk[def::I], def::I) * width(ijk[def::J], def::J) *
            width(ijk[def::K], def::K) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Width of the mesh (block) in each direction.
 *
 * \return Space_Vector (3-vector) of \e (X,Y,Z) widths in cm for the block
 * (mesh)
 */
Mesh::Space_Vector Mesh::block_widths() const
{
    REQUIRE(block_width(0) > 0.0);
    REQUIRE(block_width(1) > 0.0);
    REQUIRE(block_width(2) > 0.0);

    //creating Space_Vector to return
    return Space_Vector(block_width(0), block_width(1), block_width(2));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns the width of cell in a given direction
 *
 * \return cell width in cm
 */
double Mesh::width(size_type ijk, int dim) const
{
    REQUIRE(dim < 3);
    REQUIRE(static_cast<int>(ijk) < num_cells_dim(dim));

    ENSURE(d_width[dim][ijk] > 0);
    ENSURE(profugus::soft_equiv(d_width[dim][ijk], d_edges[dim][ijk+1]
                                 -d_edges[dim][ijk] ) );
    return d_width[dim][ijk];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns the inverse width of cell in a given direction
 *
 * \return inverse cell width in cm^{-1}
 */
double Mesh::inv_width(size_type ijk,
                       int       dim) const
{
    REQUIRE(dim < 3);
    REQUIRE(static_cast<int>(ijk) < num_cells_dim(dim));

    ENSURE(d_width[dim][ijk] > 0.0);
    return d_inv_width[dim][ijk];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the center of a cell in the specified direction.
 *
 * \param ijk \e (i,j,k) index
 * \param dim \e (I,J,K) [or \e (X,Y,Z) ] direction
 *
 * \return cell center in \e (X,Y,Z) dimension in cm
 */
double Mesh::center(size_type ijk,
                    int       dim) const
{
    REQUIRE(dim < 3);
    REQUIRE(static_cast<int>(ijk) < num_cells_dim(dim));

    return 0.5 * (d_edges[dim][ijk + 1] + d_edges[dim][ijk]);
}

} // end namespace profugus

#endif // spn_Mesh_i_hh

//---------------------------------------------------------------------------//
//              end of spn/Mesh.i.hh
//---------------------------------------------------------------------------//
