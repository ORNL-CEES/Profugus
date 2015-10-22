//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Array.i.hh
 * \author Thomas M. Evans
 * \date   Tue Dec 21 12:46:26 2010
 * \brief  Member definitions of class RTK_Array.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_RTK_Array_i_hh
#define geometry_RTK_Array_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the object in the specified position in the array.
 */
template<class T>
const typename RTK_Array<T>::Object_t& RTK_Array<T>::object(int i,
                                                            int j,
                                                            int k) const
{
    REQUIRE(id(i, j, k) >= 0 && id(i, j, k) < d_objects.size());
    REQUIRE(d_objects[id(i, j, k)]);
    return *(d_objects[id(i, j, k)]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the object in the specified position in the array.
 */
template<class T>
typename RTK_Array<T>::Object_t& RTK_Array<T>::object(int index) const
{
    REQUIRE(index >= 0 && index < d_objects.size());
    REQUIRE(d_objects[index]);
    return *(d_objects[index]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the cardinal index given logical indices in the array.
 */
template<class T>
int RTK_Array<T>::index(int i,
                        int j,
                        int k) const
{
    REQUIRE(i >= 0 && i < d_N[0]);
    REQUIRE(j >= 0 && j < d_N[1]);
    REQUIRE(k >= 0 && k < d_N[2]);
    return i + d_N[0] * (j + k * d_N[1]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the current material id.
 */
template<class T>
int RTK_Array<T>::matid(const Geo_State_t &state) const
{
    return object(state)->matid(state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the current cell id.
 */
template<class T>
int RTK_Array<T>::cellid(const Geo_State_t &state) const
{
    REQUIRE(object(state)->cellid(state) < object(state)->num_cells());
    return object(state)->cellid(state) + d_Nc_offset[
        index(state.level_coord[d_level][0],
              state.level_coord[d_level][1],
              state.level_coord[d_level][2])];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get extents in parent coordinate system.
 */
template<class T>
void RTK_Array<T>::get_extents(Space_Vector &lower, Space_Vector &upper) const
{
    using def::X; using def::Y; using def::Z;

    lower = d_corner;
    upper[X] = d_corner[X] + d_length[X];
    upper[Y] = d_corner[Y] + d_length[Y];
    upper[Z] = d_corner[Z] + d_length[Z];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transform into object coordinate system.
 */
template<class T>
typename RTK_Array<T>::Space_Vector
RTK_Array<T>::transform(const Space_Vector &r,
                        const Geo_State_t  &state) const
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(d_completed);

    REQUIRE(r[X] > d_corner[X] || soft_equiv(r[X], d_corner[X]));
    REQUIRE(r[X] < d_corner[X] + d_length[X] ||
             soft_equiv(r[X], d_corner[X] + d_length[X]));

    REQUIRE(r[Y] > d_corner[Y] || soft_equiv(r[Y], d_corner[Y]));
    REQUIRE(r[Y] < d_corner[Y] + d_length[Y] ||
            soft_equiv(r[Y], d_corner[Y] + d_length[Y]));

    REQUIRE(r[Z] > d_corner[Z] || soft_equiv(r[Z], d_corner[Z]));
    REQUIRE(r[Z] < d_corner[Z] + d_length[Z] ||
             soft_equiv(r[Z], d_corner[Z] + d_length[Z]));

#ifdef REQUIRE_ON
    // nested arrays must be initialized with lower-corners of 0.0; we only
    // allow an offset corner on the top-level array

    // Transformed coordinates and lower/upper extents
    Space_Vector lower, upper;

    // get the extents (we only need the lower)
    object(state)->get_extents(lower, upper);
    CHECK(lower[X] == 0.0);
    CHECK(lower[Y] == 0.0);
    CHECK(lower[Z] == 0.0);
#endif

    // Transformed coordinates.
    Space_Vector tr;

    // transform the coordinates to the nested array
    tr[X] = (r[X] - d_x[state.level_coord[d_level][X]]);
    tr[Y] = (r[Y] - d_y[state.level_coord[d_level][Y]]);
    tr[Z] = (r[Z] - d_z[state.level_coord[d_level][Z]]);

    return tr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get an object.
 */
template<class T>
const typename RTK_Array<T>::SP_Object&
RTK_Array<T>::object(const Geo_State_t &state) const
{
    using def::X; using def::Y; using def::Z;
    REQUIRE(d_completed);
    ENSURE(d_objects[d_layout[index(state.level_coord[d_level][X],
                                     state.level_coord[d_level][Y],
                                     state.level_coord[d_level][Z])]]);
    return d_objects[d_layout[index(state.level_coord[d_level][X],
                                    state.level_coord[d_level][Y],
                                    state.level_coord[d_level][Z])]];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update array coordinates and state when crossing a low array face.
 */
template<class T>
void RTK_Array<T>::calc_low_face(Geo_State_t &state,
                                 int          face_type,
                                 int          exiting_face)
{
    // update the coordinates in this array
    state.level_coord[d_level][face_type]--;

    // check to see if the particle has crossed this level boundary
    if (state.level_coord[d_level][face_type] < 0)
    {
        // flag indicating that the particle has crossed the low boundary of
        // this level
        state.exiting_level[d_level] = exiting_face;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update array coordinates and state when crossing a high array face.
 */
template<class T>
void RTK_Array<T>::calc_high_face(Geo_State_t &state,
                                  int          face_type,
                                  int          exiting_face)
{
    // update the coordinates in this array
    state.level_coord[d_level][face_type]++;

    // check to see if the particle has crossed this level boundary
    if (state.level_coord[d_level][face_type] > d_N[face_type] - 1)
    {
        // flag indicating that the particle has crossed the high boundary of
        // this level
        state.exiting_level[d_level] = exiting_face;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Dispatch to the correct cell counting method for objects at
 * different levels.
 */
template<class T>
int RTK_Array<T>::cell_count_dispatch(int i,
                                      int j,
                                      int k)
{
    return d_objects[id(i,j,k)]->count_cells();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add vessel to an array object.
 */
template<class T>
void RTK_Array<T>::add_vessel_to_object(int    i,
                                        int    j,
                                        int    k,
                                        double R0,
                                        double R1,
                                        double xoff,
                                        double yoff,
                                        int    vid)
{
    REQUIRE(d_objects[this->id(i, j, k)]);

    // make a (deep) copy of the object at this location
    SP_Object object(new Object_t(*d_objects[this->id(i, j, k)]));
    CHECK(object);
    CHECK(d_objects[this->id(i, j, k)] != object);

    // add this object to the list of objects
    this->id(i, j, k) = d_objects.size();
    d_objects.push_back(object);
    CHECK(d_objects[this->id(i, j, k)] == object);

    // for nested array objects, this is just a pass-through (RTK_Cells
    // require specialization)
    object->add_vessel(R0, R1, xoff, yoff, vid);
}

//---------------------------------------------------------------------------//
// RTK PIN-CELL SPECIALIZATION INLINE FUNCTIONS
//---------------------------------------------------------------------------//

template<>
inline void RTK_Array<RTK_Cell>::add_vessel_to_object(int    i,
                                                      int    j,
                                                      int    k,
                                                      double R0,
                                                      double R1,
                                                      double xoff,
                                                      double yoff,
                                                      int    vid)
{
    using def::X; using def::Y;

    REQUIRE(d_objects[this->id(i, j, k)]);

    // we have to make a new cell at the i,j,k location with the vessel
    // defined in it (otherwise we could potentially overwrite the vessel in
    // the same homogeneous cell stored at multiple locations)

    // get a reference to the existing cell
    const Object_t &cell = *d_objects[this->id(i, j, k)];
    CHECK(cell.num_regions() == 1);

    // make a new RTK_Cell that we will put at this location (but contains a
    // vessel)
    SP_Object vcell(new Object_t(cell.matid(0), cell.pitch(X), cell.pitch(Y),
                                 cell.height(), R0, R1, xoff, yoff, vid));
    CHECK(cell.num_regions() == 1);

    // add this new cell to the list of objects; the id will be the number of
    // objects currently in the vector of cells
    this->id(i, j, k) = d_objects.size();
    d_objects.push_back(vcell);
}

//---------------------------------------------------------------------------//

template<>
inline int RTK_Array<RTK_Cell>::matid(const Geo_State_t &state) const
{
    REQUIRE(d_level == 0);
    return object(state)->matid(state.region);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the current cell id.
 */
template<>
inline int RTK_Array<RTK_Cell>::cellid(const Geo_State_t &state) const
{
    REQUIRE(d_level == 0);
    REQUIRE(object(state)->cell(state.region, state.segment) <
             object(state)->num_cells());
    return object(state)->cell(state.region, state.segment) + d_Nc_offset[
        index(state.level_coord[d_level][0],
              state.level_coord[d_level][1],
              state.level_coord[d_level][2])];
}

//---------------------------------------------------------------------------//

template<>
inline RTK_Array<RTK_Cell>::Space_Vector
RTK_Array<RTK_Cell>::transform(const Space_Vector &r,
                               const Geo_State_t  &state) const
{
    using def::X; using def::Y; using def::Z;
    using std::fabs;

    REQUIRE(d_completed);
    REQUIRE(d_level == 0);

    // Transformed coordinates and lower/upper extents
    Space_Vector tr, lower, upper;

    // get the extents (we only need the lower)
    object(state)->get_extents(lower, upper);
    CHECK(lower[X] < 0.0);
    CHECK(lower[Y] < 0.0);

    // transform the coordinates to the pin cell
    tr[X] = (r[X] - d_x[state.level_coord[0][X]]) + lower[X];
    tr[Y] = (r[Y] - d_y[state.level_coord[0][Y]]) + lower[Y];
    tr[Z] = (r[Z] - d_z[state.level_coord[0][Z]]);

    ENSURE((tr[X] > lower[X] && tr[X] < upper[X]) ||
            (soft_equiv(tr[X], lower[X], 1.0e-6 * fabs(lower[X])) ||
             soft_equiv(tr[X], upper[X], 1.0e-6 * upper[X])));
    ENSURE((tr[Y] > lower[Y] && tr[Y] < upper[Y]) ||
            (soft_equiv(tr[Y], lower[Y], 1.0e-6 * fabs(lower[Y])) ||
             soft_equiv(tr[Y], upper[Y], 1.0e-6 * upper[Y])));
    return tr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Dispatch to the RTK_Cell num-cells method.
 */
template<>
inline int RTK_Array<RTK_Cell>::cell_count_dispatch(int i,
                                                    int j,
                                                    int k)
{
    return d_objects[id(i,j,k)]->num_cells();
}

//---------------------------------------------------------------------------//
// PIN-CELL SPECIALIZATION PROTOTYPES
//---------------------------------------------------------------------------//

template<>
int RTK_Array<RTK_Cell>::calc_level();

template<>
void RTK_Array<RTK_Cell>::determine_boundary_crossings(Geo_State_t &state);

template<>
void RTK_Array<RTK_Cell>::update_coordinates(const Space_Vector &r,
                                             Geo_State_t &state);

template<>
void RTK_Array<RTK_Cell>::output(std::ostream &out, int level,
                                 int obj_id) const;

} // end namespace profugus

#endif // geometry_RTK_Array_i_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_Array.i.hh
//---------------------------------------------------------------------------//
