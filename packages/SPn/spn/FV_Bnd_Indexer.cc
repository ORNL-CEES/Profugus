//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/FV_Bnd_Indexer.cc
 * \author Thomas M. Evans
 * \date   Sat Nov 24 13:46:00 2012
 * \brief  FV_Bnd_Indexer member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "FV_Bnd_Indexer.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param face face index in range [-x,+x,-y,+y,-z,+z]
 * \param N_abscissa local number of cells in the abscissa direction
 * \param N_ordinate local number of cells in the ordinate direction
 * \param G_abscissa global number of cells in the abscissa direction
 * \param G_ordinate global number of cells in the ordinate direction
 * \param bc_local local number of boundary cells on each face
 * \param bc_global global number of boundary cells one each face
 * \param off_abscissa global offset into abscissa
 * \param off_ordinate global offset into ordinate
 */
FV_Bnd_Indexer::FV_Bnd_Indexer(int face,
                               int N_abscissa,
                               int N_ordinate,
                               int bc_local[],
                               int G_abscissa,
                               int G_ordinate,
                               int bc_global[],
                               int off_abscissa,
                               int off_ordinate)
    : d_face(face)
    , d_face_off_local(0)
    , d_face_off_global(0)
{
    REQUIRE(face >= 0 && face < 6);
    REQUIRE(N_abscissa * N_ordinate == bc_local[face]);
    REQUIRE(G_abscissa * G_ordinate == bc_global[face]);

    // set local/global number of cells on the face
    d_N[0] = N_abscissa;
    d_N[1] = N_ordinate;

    d_G[0] = G_abscissa;
    d_G[1] = G_ordinate;

    // global abscissa,ordinate offsets
    d_off[0] = off_abscissa;
    d_off[1] = off_ordinate;

    // determine the local and global offsets into this face
    for (int f = 0; f < d_face; ++f)
    {
        d_face_off_local  += bc_local[f];
        d_face_off_global += bc_global[f];
    }

    ENSURE(d_face_off_local <= d_face_off_global);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of FV_Bnd_Indexer.cc
//---------------------------------------------------------------------------//
