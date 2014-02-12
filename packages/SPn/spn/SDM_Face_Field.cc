//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/SDM_Face_Field.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 30 17:05:49 2012
 * \brief  SDM_Face_Field member definitions.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "utils/Definitions.hh"
#include "SDM_Face_Field.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param mesh reference to mesh
 *
 * \param face kba::def::XYZ (X,Y,Z) or kba::def::IJK (I,J,K) enumeration
 * indicating the face type
 *
 * \sa SDM_Face_Field description.
 */
SDM_Face_Field::SDM_Face_Field(const Mesh_t &mesh,
                               int           face,
                               int           M)
    : d_face(face)
    , d_M(M)
{
    using def::X;
    using def::Y;
    using def::Z;

    Require (M > 0);

    // get the number of cells for a block face

    // for X faces: Y := abscissa and Z := ordinate
    if (face == X)
    {
        d_N_abscissa = mesh.num_cells_block_dim(Y);
        d_N_ordinate = mesh.num_cells_block_dim(Z);
    }

    // for Y faces: X := abscissa and Z := ordinate
    else if (face == Y)
    {
        d_N_abscissa = mesh.num_cells_block_dim(X);
        d_N_ordinate = mesh.num_cells_block_dim(Z);
    }

    // for Z faces: X := abscissa and Y := ordinate
    else if (face == Z)
    {
        d_N_abscissa = mesh.num_cells_block_dim(X);
        d_N_ordinate = mesh.num_cells_block_dim(Y);
    }

    else
    {
        throw profugus::assertion("Bad face type (XYZ) given.");
    }

    // total number of cells on face of block
    d_N = d_N_abscissa * d_N_ordinate;

    // size of each block matrix
    d_size = d_M * d_M;

    // allocate space
    d_field.resize(d_N * d_size);

    // fill the field
    std::fill(d_field.begin(), d_field.end(), 0.0);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                        end of SDM_Face_Field.cc
//---------------------------------------------------------------------------//
