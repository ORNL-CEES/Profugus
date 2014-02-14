//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/FV_Bnd_Indexer.hh
 * \author Thomas M. Evans
 * \date   Sat Nov 24 13:26:33 2012
 * \brief  FV_Bnd_Indexer class definition.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_FV_Bnd_Indexer_hh
#define spn_FV_Bnd_Indexer_hh

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class FV_Bnd_Indexer
 * \brief Local-to-global indexing operations for finite-volume boundary
 * conditions.
 */
/*!
 * \example spn/test/tstFV_Bnd_Indexer.cc
 *
 * Test of FV_Bnd_Indexer.
 */
//===========================================================================//

class FV_Bnd_Indexer
{
  private:
    // >>> DATA

    // Face for this boundary (-x,+x,-y,+y,-z,+z).
    int d_face;

    // Offsets into this face.
    int d_face_off_local;
    int d_face_off_global;

    // Local/global number of cells in (abscissa, ordinate).
    int d_N[2], d_G[2];

    // Global offsets into the (absicssa, ordinate).
    int d_off[2];

  public:
    // Constructor.
    FV_Bnd_Indexer(int face, int N_abscissa, int N_ordinate, int bc_local[],
                   int G_abscissa, int G_ordinate, int bc_global[],
                   int off_abscissa, int off_ordinate);

    //! Number of local cells on the face.
    int num_local() const { return d_N[0] * d_N[1]; }

    //! Number of global cells on the face.
    int num_global() const { return d_G[0] * d_G[1]; }

    //! Get local index of the boundary cell for this face.
    int local(int abscissa, int ordinate) const
    {
        Require (abscissa < d_N[0]);
        Require (ordinate < d_N[1]);
        return d_face_off_local + abscissa + ordinate * d_N[0];
    }

    //! Get the global index of the boundary cell for this face.
    int global(int abscissa, int ordinate) const
    {
        Require (abscissa < d_G[0]);
        Require (ordinate < d_G[1]);
        return d_face_off_global + abscissa + ordinate * d_G[0];
    }

    //! Local-to-global index of the boundary cell for this face.
    int l2g(int abscissa, int ordinate) const
    {
        Require (abscissa < d_N[0]);
        Require (ordinate < d_N[1]);
        return d_face_off_global + (abscissa + d_off[0]) +
            (ordinate + d_off[1]) * d_G[0];
    }
};

} // end namespace profugus

#endif // spn_FV_Bnd_Indexer_hh

//---------------------------------------------------------------------------//
//              end of spn/FV_Bnd_Indexer.hh
//---------------------------------------------------------------------------//
