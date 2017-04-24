//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Array.i.cuh
 * \author Tom Evans
 * \date   Wed Jan 04 15:43:43 2017
 * \brief  RTK_Array CUDA device class definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Array_i_cuh
#define MC_cuda_rtk_RTK_Array_i_cuh

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// DEVICE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TRACKING FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a geometry state for tracking through the array.
 */
template<class T>
__device__
void RTK_Array<T>::initialize(
    const Space_Vector &r,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    // initialize state
    state_vector.escaping_face(pid)   = Geo_State_t::NONE;
    state_vector.reflecting_face(pid) = Geo_State_t::NONE;
    state_vector.face(pid)            = Geo_State_t::NONE;
    state_vector.region(pid)          = Geo_State_t::NONE;
    state_vector.exiting_level(pid)   = {0, 0, 0};

    // find the object that this point is in
    int id = d_layout[find_object(r, state_vector, pid)];

    // Transformed coordinates.
    Space_Vector tr = transform(r, state_vector, pid);

    // initialize the state by diving into the object
    d_objects[id].initialize(tr, state_vector, pid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Track to next boundary.
 */
template<class T>
__device__
void RTK_Array<T>::distance_to_boundary(
    const Space_Vector &r,
    const Space_Vector &omega,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    // Clear reflecting face indicator
    state_vector.reflecting_face(pid) = Geo_State_t::NONE;

    // Transformed coordinates.
    Space_Vector tr = transform(r, state_vector, pid);

    // Dive through objects until we hit the pin-cell
    object(state_vector, pid).distance_to_boundary(
        tr, omega, state_vector, pid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the state of particle while in mid-track (at collision sites).
 */
template<class T>
__device__
void RTK_Array<T>::update_state(
    Geo_State_Vector_t &state_vector, int pid) const
{
    // clear reflecting face indicator
    state_vector.reflecting_face(pid) = Geo_State_t::NONE;

    // update the state of the current object if it has not escaped
    if (state_vector.escaping_face(pid) == Geo_State_t::NONE)
    {
        object(state_vector, pid).update_state(state_vector, pid);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a surface.
 */
template<class T>
__device__
void RTK_Array<T>::cross_surface(
    const Space_Vector &r,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y; using def::Z;

    // initialize exiting level flag
    state_vector.exiting_level(pid) = {0, 0, 0};

    // determine the surface crossings in each level
    determine_boundary_crossings(state_vector, pid);

    // simply return if we are at an internal pin face as everything should be
    // correctly set
    if (state_vector.exiting_face(pid) == Geo_State_t::INTERNAL)
    {
        DEVICE_CHECK(!state_vector.exiting_level(pid)[d_level]);
        return;
    }

    // update all of the downstream levels state assuming we haven't left the
    // geometry
    else if (!state_vector.exiting_level(pid)[d_level])
    {
        update_coordinates(r, state_vector, pid);
    }

    // otherwise, the particle has escaped the geometry, record the state
    else
    {
        DEVICE_CHECK(state_vector.exiting_level(pid)[d_level]);
        DEVICE_CHECK(state_vector.exiting_face(pid) > Geo_State_t::INTERNAL);

        int refl_face_index = state_vector.exiting_face(pid) -
                              Geo_State_t::MINUS_X;
        DEVICE_CHECK(refl_face_index >= 0 && refl_face_index < 6);

        // if this is a reflecting face, then reflect the particle and add 1
        // to all of the level coordinates
        if (d_reflect[refl_face_index])
        {
            // store the reflecting face
            state_vector.reflecting_face(pid) = state_vector.exiting_face(pid);

            // set the face to none representing an external pin face
            state_vector.face(pid) = Geo_State_t::NONE;

            // get correction for logical coordinates of objects in each level
            Geo_State_t::Coordinates reflect = {0, 0, 0};

            switch (state_vector.reflecting_face(pid))
            {
                case Geo_State_t::MINUS_X:
                    reflect[X] = 1;
                    break;
                case Geo_State_t::PLUS_X:
                    reflect[X] = -1;
                    break;
                case Geo_State_t::MINUS_Y:
                    reflect[Y] = 1;
                    break;
                case Geo_State_t::PLUS_Y:
                    reflect[Y] = -1;
                    break;
                case Geo_State_t::MINUS_Z:
                    reflect[Z] = 1;
                    break;
                case Geo_State_t::PLUS_Z:
                    reflect[Z] = -1;
                    break;
                default:
                    DEVICE_INSIST(false, "Lost a reflecting face.");
            }

            // loop through levels and update the coordinates to bring them
            // back into the geometry
            for (int l = 0; l <= d_level; ++l)
            {
                state_vector.level_coord(pid,l)[0] += reflect[0];
                state_vector.level_coord(pid,l)[1] += reflect[1];
                state_vector.level_coord(pid,l)[2] += reflect[2];
            }
        }

        // otherwise, the particle has left the geometry
        else
        {
            state_vector.escaping_face(pid) = state_vector.exiting_face(pid);
        }
    }
}

//---------------------------------------------------------------------------//
// ACCESSORS
//---------------------------------------------------------------------------//
/*!
 * \brief Find the object a point is in.
 */
template<class T>
__device__
int RTK_Array<T>::find_object(
    const Space_Vector &r,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y; using def::Z;

    DEVICE_REQUIRE(r[X] >= d_x.front()); DEVICE_REQUIRE(r[X] <= d_x.back());
    DEVICE_REQUIRE(r[Y] >= d_y.front()); DEVICE_REQUIRE(r[Y] <= d_y.back());
    DEVICE_REQUIRE(r[Z] >= d_z.front()); DEVICE_REQUIRE(r[Z] <= d_z.back());

    // iterators
    int                     i,j,k;
    View_Dbl::const_pointer itr, jtr, ktr;

    // find the logical indices of the object in the array
    itr = cuda::utility::lower_bound(d_x.begin(), d_x.end(), r[X]);
    i   = itr - d_x.begin() - 1;
    jtr = cuda::utility::lower_bound(d_y.begin(), d_y.end(), r[Y]);
    j   = jtr - d_y.begin() - 1;
    ktr = cuda::utility::lower_bound(d_z.begin(), d_z.end(), r[Z]);
    k   = ktr - d_z.begin() - 1;

    // check for particles on the low face of the array
    if (r[X] == d_x[0])
    {
        DEVICE_CHECK(itr == d_x.begin());

        // reset index
        i = 0;
    }
    if (r[Y] == d_y[0])
    {
        DEVICE_CHECK(jtr == d_y.begin());

        // reset index
        j = 0;
    }
    if (r[Z] == d_z[0])
    {
        DEVICE_CHECK(ktr == d_z.begin());

        // reset index
        k = 0;
    }

    DEVICE_CHECK(i >= 0 && i < d_N[X]);
    DEVICE_CHECK(j >= 0 && j < d_N[Y]);
    DEVICE_CHECK(k >= 0 && k < d_N[Z]);

    DEVICE_ENSURE(d_x[i] <= r[X] && d_x[i+1] >= r[X]);
    DEVICE_ENSURE(d_y[j] <= r[Y] && d_y[j+1] >= r[Y]);
    DEVICE_ENSURE(d_z[k] <= r[Z] && d_z[k+1] >= r[Z]);
    state_vector.level_coord(pid,d_level)[X] = i;
    state_vector.level_coord(pid,d_level)[Y] = j;
    state_vector.level_coord(pid,d_level)[Z] = k;

    return index(i, j, k);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the object on a boundary face.
 */
template<class T>
__device__
int RTK_Array<T>::find_object_on_boundary(
    const Space_Vector &r,
    int                 face,
    int                 face_type,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y; using def::Z;

    // iterators
    int                     i = 0, j = 0, k = 0;
    View_Dbl::const_pointer itr, jtr, ktr;

    // find the logical indices of the object in the array; do not search the
    // face dimension as this is known
    if (face_type != X)
    {
        itr = cuda::utility::lower_bound(d_x.begin(), d_x.end(), r[X]);
        i   = itr - d_x.begin() - 1;
        DEVICE_CHECK(d_x[i] <= r[X] && d_x[i+1] >= r[X]);
    }
    else
    {
        if (face == Geo_State_t::PLUS_X)
            i = d_N[X] - 1;
    }

    if (face_type != Y)
    {
        jtr = cuda::utility::lower_bound(d_y.begin(), d_y.end(), r[Y]);
        j   = jtr - d_y.begin() - 1;
        DEVICE_CHECK(d_y[j] <= r[Y] && d_y[j+1] >= r[Y]);
    }
    else
    {
        if (face == Geo_State_t::PLUS_Y)
            j = d_N[Y] - 1;
    }

    if (face_type != Z)
    {
        ktr = cuda::utility::lower_bound(d_z.begin(), d_z.end(), r[Z]);
        k   = ktr - d_z.begin() - 1;
        DEVICE_CHECK(d_z[k] <= r[Z] && d_z[k+1] >= r[Z]);
    }
    else
    {
        if (face == Geo_State_t::PLUS_Z)
            k = d_N[Z] - 1;
    }

    DEVICE_ENSURE(i >= 0 && i < d_N[X]);
    DEVICE_ENSURE(j >= 0 && j < d_N[Y]);
    DEVICE_ENSURE(k >= 0 && k < d_N[Z]);

    state_vector.level_coord(pid,d_level)[X] = i;
    state_vector.level_coord(pid,d_level)[Y] = j;
    state_vector.level_coord(pid,d_level)[Z] = k;

    return index(i, j, k);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the cardinal index given logical indices in the array.
 */
template<class T>
__device__
int RTK_Array<T>::index(
    int i,
    int j,
    int k) const
{
    DEVICE_REQUIRE(i >= 0 && i < d_N[0]);
    DEVICE_REQUIRE(j >= 0 && j < d_N[1]);
    DEVICE_REQUIRE(k >= 0 && k < d_N[2]);
    return i + d_N[0] * (j + k * d_N[1]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the current material id.
 */
template<class T>
__device__
int RTK_Array<T>::matid(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    return object(state_vector, pid).matid(state_vector, pid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the current cell id.
 */
template<class T>
__device__
int RTK_Array<T>::cellid(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    return object(state_vector, pid).cellid(state_vector, pid) +
           d_cell_offsets[index(state_vector.level_coord(pid,d_level)[0],
                                state_vector.level_coord(pid,d_level)[1],
                                state_vector.level_coord(pid,d_level)[2])];
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Transform into object coordinate system.
 */
template<class T>
__device__
typename RTK_Array<T>::Space_Vector
RTK_Array<T>::transform(
    const Space_Vector         &r,
    const Geo_State_Vector_t   &state_vector,
    int                         pid) const
{
    using def::X; using def::Y; using def::Z;
    using cuda::utility::soft_equiv;

    DEVICE_REQUIRE(r[X] > d_corner[X] || soft_equiv(r[X], d_corner[X]));
    DEVICE_REQUIRE(r[X] < d_corner[X] + d_length[X] ||
                   soft_equiv(r[X], d_corner[X] + d_length[X]));

    DEVICE_REQUIRE(r[Y] > d_corner[Y] || soft_equiv(r[Y], d_corner[Y]));
    DEVICE_REQUIRE(r[Y] < d_corner[Y] + d_length[Y] ||
                   soft_equiv(r[Y], d_corner[Y] + d_length[Y]));

    DEVICE_REQUIRE(r[Z] > d_corner[Z] || soft_equiv(r[Z], d_corner[Z]));
    DEVICE_REQUIRE(r[Z] < d_corner[Z] + d_length[Z] ||
                   soft_equiv(r[Z], d_corner[Z] + d_length[Z]));

    // Transformed coordinates.
    Space_Vector tr;

    // transform the coordinates to the nested array
    tr[X] = (r[X] - d_x[state_vector.level_coord(pid,d_level)[X]]);
    tr[Y] = (r[Y] - d_y[state_vector.level_coord(pid,d_level)[Y]]);
    tr[Z] = (r[Z] - d_z[state_vector.level_coord(pid,d_level)[Z]]);

    return tr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get an object.
 */
template<class T>
__device__
const typename RTK_Array<T>::Object_t&
RTK_Array<T>::object(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    using def::X; using def::Y; using def::Z;
    return d_objects[d_layout[index(state_vector.level_coord(pid,d_level)[X],
                                    state_vector.level_coord(pid,d_level)[Y],
                                    state_vector.level_coord(pid,d_level)[Z])]];
}


//---------------------------------------------------------------------------//
/*!
 * \brief Determine boundary crossings at each level starting at the lowest.
 */
template<class T>
__device__
void RTK_Array<T>::determine_boundary_crossings(
    Geo_State_Vector_t &state_vector, int pid) const
{
    using def::X; using def::Y; using def::Z;

    DEVICE_REQUIRE(d_level > 0);

    // dive into the object and see if it crosses a boundary
    object(state_vector, pid).determine_boundary_crossings(state_vector, pid);

    // process particles that cross an array boundary on the previous level,
    // they may cross a boundary at this level as well
    if (state_vector.exiting_level[d_level - 1])
    {
        switch (state_vector.exiting_face(pid))
        {
            case Geo_State_t::MINUS_X:
                calc_low_face(state_vector, pid, X, Geo_State_t::MINUS_X);
                break;
            case Geo_State_t::PLUS_X:
                calc_high_face(state_vector, pid, X, Geo_State_t::PLUS_X);
                break;
            case Geo_State_t::MINUS_Y:
                calc_low_face(state_vector, pid, Y, Geo_State_t::MINUS_Y);
                break;
            case Geo_State_t::PLUS_Y:
                calc_high_face(state_vector, pid, Y, Geo_State_t::PLUS_Y);
                break;
            case Geo_State_t::MINUS_Z:
                calc_low_face(state_vector, pid, Z, Geo_State_t::MINUS_Z);
                break;
            case Geo_State_t::PLUS_Z:
                calc_high_face(state_vector, pid, Z, Geo_State_t::PLUS_Z);
                break;
            default:
                DEVICE_INSIST(
                    false, "Not at a valid array surface crossing.");
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the array coordinates in all levels where the particle has
 * crossed a face.
 */
template<class T>
__device__
void RTK_Array<T>::update_coordinates(
    const Space_Vector &r,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y; using def::Z;

    DEVICE_REQUIRE(d_level > 0);

    // transform coordinates into object's coordinate system
    Space_Vector tr = transform(r, state_vector, pid);

    // if the child object has exited the level, then update the coordinates
    if (state_vector.exiting_level(pid)[d_level - 1])
    {
        switch (state_vector.exiting_face(pid))
        {
            // update the coordinates of the object across the given face
            case Geo_State_t::MINUS_X:
                object(state_vector, pid).find_object_on_boundary(
                    tr, Geo_State_t::PLUS_X, X, state_vector, pid);
                break;
            case Geo_State_t::PLUS_X:
                object(state_vector, pid).find_object_on_boundary(
                    tr, Geo_State_t::MINUS_X, X, state_vector, pid);
                break;
            case Geo_State_t::MINUS_Y:
                object(state_vector, pid).find_object_on_boundary(
                    tr, Geo_State_t::PLUS_Y, Y, state_vector, pid);
                break;
            case Geo_State_t::PLUS_Y:
                object(state_vector, pid).find_object_on_boundary(
                    tr, Geo_State_t::MINUS_Y, Y, state_vector, pid);
                break;
            case Geo_State_t::MINUS_Z:
                object(state_vector, pid).find_object_on_boundary(
                    tr, Geo_State_t::PLUS_Z, Z, state_vector, pid);
                break;
            case Geo_State_t::PLUS_Z:
                object(state_vector, pid).find_object_on_boundary(
                    tr, Geo_State_t::MINUS_Z, Z, state_vector, pid);
                break;
        }
    }

    // go to child object and recursively update all coordinates
    object(state_vector, pid).update_coordinates(tr, state_vector, pid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update array coordinates and state when crossing a low array face.
 */
template<class T>
__device__
void RTK_Array<T>::calc_low_face(
    Geo_State_Vector_t &state_vector,
    int                 pid
    int                 face_type,
    int                 exiting_face) const
{
    // update the coordinates in this array
    --state_vector.level_coord(pid,d_level)[face_type];

    // check to see if the particle has crossed this level boundary
    if (state_vector.level_coord(pid,d_level)[face_type] < 0)
    {
        // flag indicating that the particle has crossed the low boundary of
        // this level
        state_vector.exiting_level(pid)[d_level] = exiting_face;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update array coordinates and state when crossing a high array face.
 */
template<class T>
__device__
void RTK_Array<T>::calc_high_face(
    Geo_State_Vector_t &state_vector,
    int                 pid
    int                 face_type,
    int                 exiting_face) const
{
    // update the coordinates in this array
    ++state_vector.level_coord(pid,d_level)[face_type];

    // check to see if the particle has crossed this level boundary
    if (state_vector.level_coord(pid,d_level)[face_type] > d_N[face_type] - 1)
    {
        // flag indicating that the particle has crossed the high boundary of
        // this level
        state_vector.exiting_level(pid)[d_level] = exiting_face;
    }
}

//---------------------------------------------------------------------------//
// SPECIALIZATIONS ON RTK_CELL
//---------------------------------------------------------------------------//

template<>
__device__
inline RTK_Array<RTK_Cell>::Space_Vector
RTK_Array<RTK_Cell>::transform(
    const Space_Vector        &r,
    const Geo_State_Vector_t  &state_vector,
    int                        pid) const
{
    using def::X; using def::Y; using def::Z;
    using cuda::utility::soft_equiv;
    using std::fabs;

    DEVICE_REQUIRE(d_level == 0);

    // Transformed coordinates and lower/upper extents
    Space_Vector tr, lower, upper;

    // get the extents (we only need the lower)
    object(state_vector, pid).get_extents(lower, upper);
    DEVICE_CHECK(lower[X] < 0.0);
    DEVICE_CHECK(lower[Y] < 0.0);

    // transform the coordinates to the pin cell
    tr[X] = (r[X] - d_x[state_vector.level_coord(pid,0)[X]]) + lower[X];
    tr[Y] = (r[Y] - d_y[state_vector.level_coord(pid,0)[Y]]) + lower[Y];
    tr[Z] = (r[Z] - d_z[state_vector.level_coord(pid,0)[Z]]);

    DEVICE_ENSURE((tr[X] > lower[X] && tr[X] < upper[X]) ||
                  (soft_equiv(tr[X], lower[X], 1.0e-6 * fabs(lower[X])) ||
                   soft_equiv(tr[X], upper[X], 1.0e-6 * upper[X])));
    DEVICE_ENSURE((tr[Y] > lower[Y] && tr[Y] < upper[Y]) ||
                  (soft_equiv(tr[Y], lower[Y], 1.0e-6 * fabs(lower[Y])) ||
                   soft_equiv(tr[Y], upper[Y], 1.0e-6 * upper[Y])));
    return tr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the current cell id.
 */
template<>
__device__
inline int RTK_Array<RTK_Cell>::cellid(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    DEVICE_REQUIRE(d_level == 0);
    return object(state_vector, pid).cell(state_vector.region(pid),
                                          state_vector.segment(pid)) +
            d_cell_offsets[index(state_vector.level_coord(pid,d_level)[0],
                                 state_vector.level_coord(pid,d_level)[1],
                                 state_vector.level_coord(pid,d_level)[2])];
}

//---------------------------------------------------------------------------//

template<>
__device__
inline int RTK_Array<RTK_Cell>::matid(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    DEVICE_REQUIRE(d_level == 0);
    return object(state_vector, pid).matid(state_vector.region(pid));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine boundary crossing in array of pin cells.
 */
template<>
__device__
inline void RTK_Array<RTK_Cell>::determine_boundary_crossings(
    Geo_State_Vector_t &state_vector, int pid) const
{
    using def::X; using def::Y; using def::Z;

    DEVICE_REQUIRE(d_level == 0);

    switch (state_vector.exiting_face(pid))
    {
        // process internal pin-cell surface crossings
        case Geo_State_t::INTERNAL:
            // call the pin-cell surface crossing routine for internal surfaces
            object(state_vector, pid).cross_surface(state_vector, pid);
            break;
        // otherwise process particles leaving the pin cell by updating the
        // region id to the moderator region in the next pin cell
        case Geo_State_t::MINUS_X:
            calc_low_face(state_vector, pid, X, Geo_State_t::MINUS_X);
            break;
        case Geo_State_t::PLUS_X:
            calc_high_face(state_vector, pid, X, Geo_State_t::PLUS_X);
            break;
        case Geo_State_t::MINUS_Y:
            calc_low_face(state_vector, pid, Y, Geo_State_t::MINUS_Y);
            break;
        case Geo_State_t::PLUS_Y:
            calc_high_face(state_vector, pid, Y, Geo_State_t::PLUS_Y);
            break;
        case Geo_State_t::MINUS_Z:
            calc_low_face(state_vector, pid, Z, Geo_State_t::MINUS_Z);
            break;
        case Geo_State_t::PLUS_Z:
            calc_high_face(state_vector, pid, Z, Geo_State_t::PLUS_Z);
            break;
        default:
            // otherwise we aren't at a surface crossing at all!
            DEVICE_INSIST(false, "Not at a valid pin-cell surface crossing.");
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the array coordinates in a pin-cell array if the particle
 * leaves the pin cell.
 */
template<>
inline void RTK_Array<RTK_Cell>::update_coordinates(
    const Space_Vector &r,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    DEVICE_REQUIRE(d_level == 0);

    // if we exited the last pin-cell array, set the region coordinates in the
    // new array; for radial entrance (x,y), we must be entering the moderator
    // region
    if (state_vector(pid).exiting_face != Geo_State_t::INTERNAL)
    {
        DEVICE_CHECK(state_vector(pid).exiting_face != Geo_State_t::NONE);
        DEVICE_CHECK(state_vector.region(pid) == Geo_State_t::VESSEL ?
                     object(state_vector, pid).has_vessel() : true);

        // update the face of the state
        state_vector.face(pid) = Geo_State_t::NONE;

        // set to moderator in the pin cell
        state_vector.region(pid) = object(state_vector, pid).num_regions() - 1;

        // need to transport into coordinate system of pin cell
        Space_Vector tr = transform(r, state_vector, pid);

        // calculate the segment in the new pin-cell
        state_vector.segment(pid) =
            object(state_vector, pid).segment(tr[0], tr[1]);
        DEVICE_CHECK(state_vector.segment(pid) <
                     object(state_vector, pid).num_segments());

        // update if crossing into pin cell from high or low Z-face or if the
        // adjoining cell has a vessel
        if (state_vector.exiting_face(pid) == Geo_State_t::MINUS_Z ||
            state_vector.exiting_face(pid) == Geo_State_t::PLUS_Z  ||
            object(state_vector, pid).has_vessel())
        {
            // determine the region of the point entering the pin-cell through
            // the Z-face
            state_vector.region(pid) =
                object(state_vector, pid).region(tr[0], tr[1]);
        }
    }
}

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_Array_i_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.i.cuh
//---------------------------------------------------------------------------//
