//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_geometry/Mesh_State.hh
 * \author Steven Hamilton
 * \date   Monday July 21 18:16:58 2014
 * \brief  Mesh_State class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_geometry_Mesh_State_hh
#define cuda_geometry_Mesh_State_hh

#include "utils/Definitions.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \struct Mesh_State
 * \brief Track particle location on a cartesian mesh
 */
//===========================================================================//


struct Mesh_State
{
  public:

      typedef cuda::Coordinates  Coordinates;
      typedef cuda::Space_Vector Space_Vector;

    //! Faces in pin-cell.
    enum Face {
        MINUS_X = 0,
        PLUS_X ,
        MINUS_Y,
        PLUS_Y ,
        MINUS_Z,
        PLUS_Z,
        END_FACES
    };

  public:
    // >>> PUBLIC DATA

    //! Indices along the mesh grid if inside (invalid if not)
    Coordinates ijk;

    //! Position
    Space_Vector d_r;

    //! Direction
    Space_Vector d_dir;

    //! Coordinates of next cell (not pickled)
    Coordinates next_ijk;

    //! Distance to next cell
    double next_dist;
};

} // end namespace profugus

#endif // cuda_geometry_Mesh_State_hh

//---------------------------------------------------------------------------//
//              end of cuda_geometry/Mesh_State.hh
//---------------------------------------------------------------------------//
