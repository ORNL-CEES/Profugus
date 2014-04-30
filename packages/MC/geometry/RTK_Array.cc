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

    Require (d_level == 0);

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
            Validate(false, "Not at a valid pin-cell surface crossing.");
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
    Require (d_level == 0);

    // if we exited the last pin-cell array, set the region coordinates in the
    // new array; for radial entrance (x,y), we must be entering the moderator
    // region
    if (state.exiting_face != Geo_State_t::INTERNAL)
    {
        Check (state.exiting_face != Geo_State_t::NONE);
        Check (state.region == Geo_State_t::VESSEL ?
               object(state)->has_vessel() : true);

        // update the face of the state
        state.face = Geo_State_t::NONE;

        // set to moderator in the pin cell
        state.region = object(state)->num_regions() - 1;

        // need to transport into coordinate system of pin cell
        Space_Vector tr = transform(r, state);

        // calculate the segment in the new pin-cell
        state.segment = object(state)->segment(tr[0], tr[1]);
        Check (state.segment < object(state)->num_segments());

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
 * \brief Set symmetric/mapped cell ids
 *
 * \param this_offset this object's cell offset
 * \param map_offset map object's cell offset
 * \param t_it Iterator into d_mapped_cells for this object
 * \param m_it Iterator into d_mapped_cells for mapped object
 *
 */
template<>
void RTK_Array<RTK_Cell>::set_mapped_cells(int               this_offset,
                                           int               map_offset,
                                           Vec_Int::iterator t_it,
                                           Vec_Int::iterator m_it)
{
    Check( d_level == 0 );

    int n_blocks = d_N[def::X] * d_N[def::Y];

    Vec_Int mapped_array(n_blocks);

    int index = 0;
    int mapped_index;
    for (int j = 0; j < d_N[def::Y]; ++j)
    {
        mapped_index = d_N[def::X] * d_N[def::Y] - 1 - j;
        for (int i = 0; i < d_N[def::X]; ++i, ++index)
        {
            // setting mapped array
            mapped_array[index] = mapped_index;

            // Checking mapped_index and index have same number of cells
            Check( d_objects[ d_layout[index       ]]->num_cells() ==
                   d_objects[ d_layout[mapped_index]]->num_cells() );

            // decrementing mapped_index
            mapped_index -= d_N[def::X];
        }
    }

    // local cell index and local mapped cell index
    int local_index, map_index;
    // number of segments and rings/shells for each object
    int num_segments, num_rings;

    for (int i = 0; i < n_blocks; ++i)
    {
        // num segments
        num_segments = d_objects[ d_layout[i] ]->num_segments();
        num_rings    = d_objects[ d_layout[i] ]->num_shells();

        // segment offset
        // ( 0 and 3 are switched) (1 and 2 are not)
        Vec_Int offs( num_segments, 0 );
        if ( num_segments == 4)
        {
            offs[0] =  3 * ( num_rings + 1);
            offs[3] = -3 * ( num_rings + 1);
        }

        // 1) "Block" is completely independent (i.e mapped_index != index)
        // 2) "Block" is split (i.e mapped_index == index )
        // Assumes only up to 4 segments

        // In case 2), we repeat the assigment for segments 1 and 2 in the
        // four segment case because m_off == t_off and local_index =
        // map_index.  In segments 0 and 3 local_index != map_index

        // loop over segments
        for (int s = 0; s < num_segments; ++s)
        {
            // loop over rings
            for (int r = 0; r < num_rings + 1; ++r)
            {
                // local index
                local_index = s * (num_rings + 1) + r;

                // mapped cell's local index
                map_index = local_index + offs[s];

                // set this cell
                *(t_it + local_index + d_Nc_offset[i])
                    = map_offset + map_index + d_Nc_offset[ mapped_array[i] ];

                // set mapped cell
                *(m_it + map_index + d_Nc_offset[ mapped_array[i] ] )
                    = this_offset + local_index + d_Nc_offset[ i ];

            }
        }
    } // end loop over objects

} // end set_mapped_cells(int, int, Vec_Int::iterator, Vec_Int::iterator)

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

    Require (d_level == 0);

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
