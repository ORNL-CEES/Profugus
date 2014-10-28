//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/geometry/RTK_State.hh
 * \author Thomas M. Evans
 * \date   Tuesday April 29 16:9:53 2014
 * \brief  RTK_State struct definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_geometry_RTK_State_hh
#define core_geometry_RTK_State_hh

#include <Utils/config.h>
#include "utils/Vector_Lite.hh"
#include "utils/Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \struct RTK_State
 * \brief  Handle to basic Reactor ToolKit pin-cell, core geometry package.
 *
 * The RTK_State is a handle into the basic RTK core geometry package that
 * describes the position and state of a particle at any point in time.
 */
//===========================================================================//

struct RTK_State
{
    //! Position.
    def::Space_Vector d_r;

    //! Direction.
    def::Space_Vector d_dir;

    // >>> REQUIRED DEFINITIONS

    // Pack the geometric state.
    static int packed_bytes() { return (16 * SIZEOF_INT + 6 * SIZEOF_DOUBLE); }
    void pack(char *buffer) const;

    // Unpack the geometric state.
    void unpack(const char *buffer);

    // >>> GEOMETRIC STATE

    //! Faces in pin-cell and vessel.
    enum Faces {MODERATOR = 0,
                NONE      = 900,
                INTERNAL  = 901,
                MINUS_X   = 1000,
                PLUS_X    = 1001,
                MINUS_Y   = 1002,
                PLUS_Y    = 1003,
                MINUS_Z   = 1004,
                PLUS_Z    = 1005,
                R0_VESSEL = 2000,
                R1_VESSEL = 2001,
                VESSEL    = 2002};
    static const int plus_face[3];
    static const int minus_face[3];

    //@{
    //! Pin-cell semantics.
    int    region;
    int    segment;
    int    face;
    int    next_region;
    int    next_segment;
    int    next_face;
    double dist_to_next_region;
    //@}

    //! Exiting face indicator.
    int exiting_face;

    //! Max levels supported.
    static const int max_levels = 3;

    //! Coordinates in array at each level.
    Vector_Lite<Vector_Lite<int, 3>, max_levels> level_coord;

    //! Crossing boundary indicator by level.
    Vector_Lite<int, max_levels> exiting_level;

    //! Escaping face in geometry.
    int escaping_face;

    //! Reflecting face in geometry.
    int reflecting_face;
};

} // end namespace profugus

#endif // core_geometry_RTK_State_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_State.hh
//---------------------------------------------------------------------------//
