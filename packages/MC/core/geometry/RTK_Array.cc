//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Array.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 21 12:46:26 2010
 * \brief  RTK_Array specializations on RTK_Cell.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iomanip>
#include <algorithm>

#include "harness/Soft_Equivalence.hh"
#include "RTK_Array.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// RTK PIN-CELL SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Specialization on RTK_Cell.
 */
template<>
int RTK_Array<RTK_Cell>::calc_level()
{
    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine boundary crossing in array of pin cells.
 */
template<>
void RTK_Array<RTK_Cell>::determine_boundary_crossings(Geo_State_t &state)
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(d_level == 0);

    switch (state.exiting_face)
    {
        // process internal pin-cell surface crossings
        case Geo_State_t::INTERNAL:
            // call the pin-cell surface crossing routine for internal surfaces
            object(state)->cross_surface(state);
            break;
        // otherwise process particles leaving the pin cell by updating the
        // region id to the moderator region in the next pin cell
        case Geo_State_t::MINUS_X:
            calc_low_face(state, X, Geo_State_t::MINUS_X);
            break;
        case Geo_State_t::PLUS_X:
            calc_high_face(state, X, Geo_State_t::PLUS_X);
            break;
        case Geo_State_t::MINUS_Y:
            calc_low_face(state, Y, Geo_State_t::MINUS_Y);
            break;
        case Geo_State_t::PLUS_Y:
            calc_high_face(state, Y, Geo_State_t::PLUS_Y);
            break;
        case Geo_State_t::MINUS_Z:
            calc_low_face(state, Z, Geo_State_t::MINUS_Z);
            break;
        case Geo_State_t::PLUS_Z:
            calc_high_face(state, Z, Geo_State_t::PLUS_Z);
            break;
        default:
            // otherwise we aren't at a surface crossing at all!
            VALIDATE(false, "Not at a valid pin-cell surface crossing.");
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the array coordinates in a pin-cell array if the particle
 * leaves the pin cell.
 */
template<>
void RTK_Array<RTK_Cell>::update_coordinates(const Space_Vector &r,
                                             Geo_State_t        &state)
{
    REQUIRE(d_level == 0);

    // if we exited the last pin-cell array, set the region coordinates in the
    // new array; for radial entrance (x,y), we must be entering the moderator
    // region
    if (state.exiting_face != Geo_State_t::INTERNAL)
    {
        CHECK(state.exiting_face != Geo_State_t::NONE);
        CHECK(state.region == Geo_State_t::VESSEL ?
               object(state)->has_vessel() : true);

        // update the face of the state
        state.face = Geo_State_t::NONE;

        // set to moderator in the pin cell
        state.region = object(state)->num_regions() - 1;

        // need to transport into coordinate system of pin cell
        Space_Vector tr = transform(r, state);

        // calculate the segment in the new pin-cell
        state.segment = object(state)->segment(tr[0], tr[1]);
        CHECK(state.segment < object(state)->num_segments());

        // update if crossing into pin cell from high or low Z-face or if the
        // adjoining cell has a vessel
        if (state.exiting_face == Geo_State_t::MINUS_Z ||
            state.exiting_face == Geo_State_t::PLUS_Z  ||
            object(state)->has_vessel())
        {
            // determine the region of the point entering the pin-cell through
            // the Z-face
            state.region = object(state)->region(tr[0], tr[1]);
        }
    }
} //end update_coordinates( const Space_Vector, Geo_State_t )

//---------------------------------------------------------------------------//
/*!
 * \brief Output implementation pincells arrays
 */
template<>
void RTK_Array<RTK_Cell>::output(std::ostream &out,
                                 int           level,
                                 int           obj_id) const
{
    using std::endl; using std::setw; using std::setfill; using std::left;
    using std::right; using std::fixed; using std::setprecision;
    using def::I; using def::J; using def::K;

    REQUIRE(d_level == 0);

    out << "----------------------------------------"
        << "----------------------------------------" << endl;
    out << endl;
    out << left;

    out << "Output pin-cell array " << setw(3) << obj_id << " from level "
        << setw(3) << level << endl;
    out << endl;

    for (int k = 0; k < d_N[K]; ++k)
    {
        out << "- Array at k = " << setw(4) << k << fixed
            << setprecision(4) << setw(12) << d_z[k] << endl;
        for (int j = d_N[J]-1; j > -1; --j)
        {
            for (int i = 0; i < d_N[I]; ++i)
            {
                out << setw(5) << id(i, j, k);
            }
            out << endl;
        }
    }
    out << endl;

    // output nested objects
    for (int n = 0; n < d_objects.size(); ++n)
    {
        if (d_objects[n])
        {
            d_objects[n]->output(out, obj_id, n);
        }
    }
    out << endl;
    out << right;
    out << "----------------------------------------"
        << "----------------------------------------" << endl;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RTK_Array.cc
//---------------------------------------------------------------------------//
