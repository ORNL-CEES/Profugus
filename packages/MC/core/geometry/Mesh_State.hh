//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/geometry/Mesh_State.hh
 * \author Thomas M. Evans and Seth R Johnson
 * \date   Monday July 21 18:16:58 2014
 * \brief  Mesh_State class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_geometry_Mesh_State_hh
#define core_geometry_Mesh_State_hh

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Vector_Lite.hh"

namespace profugus
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
    typedef profugus::Vector_Lite<int, 3> Coordinates;
    typedef def::Space_Vector             Space_Vector;

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
    def::Space_Vector pos;

    //! Direction
    def::Space_Vector dir;

    //! Coordinates of next cell (not pickled)
    Coordinates next_ijk;

    //! Distance to next cell
    double next_dist;

  public:
    // >>> PACKING
    // Pack the geometric state.
    static int packed_bytes() { return (4 * SIZEOF_INT + 6 * SIZEOF_DOUBLE); }
    void pack(char *buffer) const;

    // Unpack the geometric state.
    void unpack(const char *buffer);
};

} // end namespace profugus

#endif // core_geometry_Mesh_State_hh

//---------------------------------------------------------------------------//
//              end of geometry/Mesh_State.hh
//---------------------------------------------------------------------------//
