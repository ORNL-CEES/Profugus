//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_State.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 25 14:38:37 2011
 * \brief  RTK_State member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "utils/Packing_Utils.hh"
#include "RTK_State.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// STATIC MEMBERS
//---------------------------------------------------------------------------//

const int RTK_State::plus_face[3]  = {PLUS_X, PLUS_Y, PLUS_Z};
const int RTK_State::minus_face[3] = {MINUS_X, MINUS_Y, MINUS_Z};

//---------------------------------------------------------------------------//
// PACK/UNPACK FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \return Pack the state into a char buffer.
 *
 * \param buffer char buffer that should have at least size
 * RTK_State::packed_bytes() available.
 */
void RTK_State::pack(char *buffer) const
{
    Require (buffer);
    Require (escaping_face == NONE);

    // make a packer
    Packer p;
    p.set_buffer(packed_bytes(), buffer);
    Check (sizeof(level_coord) == max_levels * 3 * SIZEOF_INT);

    // pack the data
    p << d_r << d_dir << region << segment << face << next_region
      << next_segment << next_face << exiting_face << level_coord;

    Ensure (p.get_ptr() == p.end());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack a state.
 *
 * \param buffer char buffer with a packed state (from pack()) that has size
 * packed_bytes()
 */
void RTK_State::unpack(const char *buffer)
{
    Require (buffer);

    // make an unpacker
    Unpacker u;
    u.set_buffer(packed_bytes(), buffer);

    // unpack the data
    u >> d_r >> d_dir >> region >> segment >> face >> next_region
      >> next_segment >> next_face >> exiting_face >> level_coord;

    // initialize the escaping face to none as an escaped particle should
    // never be packed
    escaping_face = NONE;

    Ensure (u.get_ptr() == u.end());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RTK_State.cc
//---------------------------------------------------------------------------//
