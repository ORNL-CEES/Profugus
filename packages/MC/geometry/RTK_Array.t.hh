//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Array.t.hh
 * \author Thomas M. Evans
 * \date   Tue Dec 21 12:46:26 2010
 * \brief  RTK_Array template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_RTK_Array_t_hh
#define geometry_RTK_Array_t_hh

#include <iomanip>

#include "harness/Soft_Equivalence.hh"
#include "utils/Definitions.hh"
#include "RTK_Array.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
template<class T>
RTK_Array<T>::RTK_Array(int Nx,
                        int Ny,
                        int Nz,
                        int num_objects)
    : d_N(Nx, Ny, Nz)
    , d_layout(d_N[0] * d_N[1] * d_N[2], 0)
    , d_objects(num_objects)
    , d_x(Nx + 1, 0.0)
    , d_y(Ny + 1, 0.0)
    , d_z(Nz + 1, 0.0)
    , d_reflect(6, 0)
    , d_num_cells(d_N[0] * d_N[1] * d_N[2], 0)
    , d_Nc_offset(d_N[0] * d_N[1] * d_N[2] + 1, 0)
    , d_total_cells(0)
    , d_vessel(false)
    , d_completed(false)
{
    using def::X; using def::Y; using def::Z;

    // calculate the level for quick access
    d_level = calc_level();

    ENSURE(d_level < Geo_State_t::max_levels);
    ENSURE(d_layout.size() == size());
    ENSURE(size() > 0);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Assign object to an array id.
 */
template<class T>
void RTK_Array<T>::assign_object(SP_Object object,
                                 int       id)
{
    REQUIRE(object);
    REQUIRE(id >= 0 && id < d_objects.size());
    d_objects[id] = object;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set reflecting faces.
 *
 * \param reflecting_faces a 6-dimension vector ordered \c (x-,x+,y-,y+,z-,z+)
 * where a zero indicates a vacuum (default) boundary and 1 indicates a
 * reflecting boundary
 *
 * \pre !complete()
 */
template<class T>
void RTK_Array<T>::set_reflecting(const Vec_Int &reflecting_faces)
{
    REQUIRE(!d_completed);
    REQUIRE(reflecting_faces.size() == 6);

    // assign the reflecting faces
    d_reflect = reflecting_faces;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set an outer vessel.
 *
 * \param r0 inner radius of vessel
 * \param r1 outer radius of vessel
 * \param id material id of vessel
 *
 * The vessel origin is at the geometric centerpoint of the array.
 *
 * \pre (r0,r1) < (x,y) and !complete()
 */
template<class T>
void RTK_Array<T>::set_vessel(double r0,
                              double r1,
                              int    vid)
{
    REQUIRE(!d_completed);

    // set vessel on to true
    d_vessel = true;

    // set the radii
    d_r[0] = r0;
    d_r[1] = r1;

    // set the vessel material id
    d_vessel_id = vid;

    ENSURE(d_r[1] > d_r[0]);
    ENSURE(d_vessel_id >= 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Complete the construction of the array.
 *
 * \param low_x low x coordinate of the array
 * \param low_y low y coordinate of the array
 * \param low_z low z coordinate of the array
 */
template<class T>
void RTK_Array<T>::complete(double low_x,
                            double low_y,
                            double low_z)
{
    using def::X; using def::Y; using def::Z;

    // assign the corner point
    d_corner = Space_Vector(low_x, low_y, low_z);

    // array id
    int id = 0;

    // loop through ids and figure out dimensions

    // X-dimensions (look at (j,k) = (0,0) row)
    d_x[0] = low_x;
    for (int i = 0; i < d_N[X]; ++i)
    {
        // array id
        id = d_layout[index(i, 0, 0)];
        CHECK(id >= 0 && id < d_objects.size());
        CHECK(d_objects[id]);
        CHECK(d_objects[id]->completed());

        // assign the lengths
        d_x[i+1] = d_x[i] + d_objects[id]->pitch(X);
    }

    // Y-dimensions (look at (i,k) = (0,0) row)
    d_y[0] = low_y;
    for (int j = 0; j < d_N[Y]; ++j)
    {
        // array id
        id = d_layout[index(0, j, 0)];
        CHECK(id >= 0 && id < d_objects.size());
        CHECK(d_objects[id]);
        CHECK(d_objects[id]->completed());

        // assign the lengths
        d_y[j+1] = d_y[j] + d_objects[id]->pitch(Y);
    }

    // Z-dimensions (look at (i,j) = (0,0) row)
    d_z[0] = low_z;
    for (int k = 0; k < d_N[Z]; ++k)
    {
        // array id
        id = d_layout[index(0, 0, k)];
        CHECK(id >= 0 && id < d_objects.size());
        CHECK(d_objects[id]);
        CHECK(d_objects[id]->completed());

        // assign the lengths
        d_z[k+1] = d_z[k] + d_objects[id]->height();
    }

    // calculate lengths
    d_length[X] = d_x.back() - low_x;
    d_length[Y] = d_y.back() - low_y;
    d_length[Z] = d_z.back() - low_z;

    // count cells at each level
    count_cells();

    // make offsets for each array element giving the first cell index in the
    // object
    std::transform(d_num_cells.begin(), d_num_cells.end(), d_Nc_offset.begin(),
                   d_Nc_offset.begin() + 1, std::plus<int>());
    CHECK(d_Nc_offset.back() == d_total_cells);

    // setup core vessel if activated
    if (d_vessel)
    {
        CHECK(2.0 * d_r[1] <= d_length[X]);
        CHECK(2.0 * d_r[1] <= d_length[Y]);

        // determine the origin
        d_origin[X] = (d_x.back() + d_x.front()) * 0.5;
        d_origin[Y] = (d_y.back() + d_y.front()) * 0.5;

        VALIDATE(d_origin[X] + d_r[1] <= d_x.back() &&
                  d_origin[Y] + d_r[1] <= d_y.back(),
                  "Radius of core vessel to large: "
                  << d_origin[X] + d_r[1] << " > " << d_x.back()
                  << " or " << d_origin[Y] + d_r[1] << " > " << d_y.back());

        // offsets (from origin of vessel to left edge of object)
        double xoff = 0.0, yoff = 0.0;

        // vector of object coordinates that have a vessel
        Vec_Int_Pair vobj;

        // near/far radii that bound the cell
        double near = 0.0, far = 0.0;

        // the vessel radii squared
        double r02 = d_r[0]*d_r[0];
        double r12 = d_r[1]*d_r[1];

        // add vessel to all objects; start by investigating objects in the PP
        // quadrant
        int start_x = d_N[X] / 2, start_y = d_N[Y] / 2;
        for (int j = start_y; j < d_N[Y]; ++j)
        {
            for (int i = start_x; i < d_N[X]; ++i)
            {
                // calculate the near/far radii that bound the cell (radii
                // squared)
                near = (d_x[i]-d_origin[X]) * (d_x[i]-d_origin[X]) +
                       (d_y[j]-d_origin[Y]) * (d_y[j]-d_origin[Y]);
                far  = (d_x[i+1]-d_origin[X]) * (d_x[i+1]-d_origin[X]) +
                       (d_y[j+1]-d_origin[Y]) * (d_y[j+1]-d_origin[Y]);

                // check to see if the vessel bisects the cell
                if (r02 > far || r12 < near)
                {
                    continue;
                }

                // now check to see if one of the vessel radii bisects the
                // cell
                if (r02 > near || r12 < far)
                {
                    // add this object to the list of array objects that are
                    // bisected by the vessel
                    vobj.push_back(std::pair<int,int>(i, j));
                }
            }
        }

        // add the symmetric parts
        Vec_Int_Pair scells;

        // update starting index to account for even/odd number of objects
        start_x += d_N[X] % 2;
        start_y += d_N[Y] % 2;

        // flip to top-left quadrant
        for (Vec_Int_Pair_Itr p = vobj.begin(); p != vobj.end(); ++p)
        {
            // add symmetric cells in the top-left quadrant
            if (p->first >= start_x)
            {
                scells.push_back(std::pair<int,int>(
                                     d_N[X] - p->first - 1, p->second));
            }
        }

        // add the symmetric cells to the list
        vobj.insert(vobj.end(), scells.begin(), scells.end());
        scells.clear();

        // flip to the lower half
        for (Vec_Int_Pair_Itr p = vobj.begin(); p != vobj.end(); ++p)
        {
            // add symmetric cells in the bottom half
            if (p->second >= start_y)
            {
                scells.push_back(std::pair<int,int>(
                                     p->first, d_N[Y] - p->second - 1));
            }
        }

        // add the symmetric cells to the list
        vobj.insert(vobj.end(), scells.begin(), scells.end());;
        scells.clear();

        // recursively drift through each object and add the vessel
        for (Vec_Int_Pair_Itr p = vobj.begin(); p != vobj.end(); ++p)
        {
            CHECK(d_objects[this->id(p->first, p->second, 0)]);

            // determine offsets to this object (the offsets can be negative;
            // the are defined relative to the lower-left corner of each
            // object array)
            xoff = d_x[p->first]  - d_origin[X];
            yoff = d_y[p->second] - d_origin[Y];

            // add the vessel to the nested object
            for (int k = 0; k < d_N[Z]; ++k)
            {
                add_vessel_to_object(p->first, p->second, k, d_r[0], d_r[1],
                                     xoff, yoff, d_vessel_id);
            }
        }
    }

    // set complete to true
    d_completed = true;

#ifdef REMEMBER_ON
    for (int k = 0; k < d_N[Z]; ++k)
    {
        for (int j = 0; j < d_N[Y]; ++j)
        {
            for (int i = 0; i < d_N[X]; ++i)
            {
                // array id
                id = d_layout[index(i, j, k)];

                // check the lengths
                ENSURE(soft_equiv(d_objects[id]->pitch(X), dx(i)));
                ENSURE(soft_equiv(d_objects[id]->pitch(Y), dy(j)));
                ENSURE(soft_equiv(d_objects[id]->height(), dz(k)));
            }
        }
    }
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a geometry state for tracking through the array.
 */
template<class T>
void RTK_Array<T>::initialize(const Space_Vector &r,
                              Geo_State_t        &state) const
{
    REQUIRE(d_completed);

    // initialize state
    state.escaping_face   = Geo_State_t::NONE;
    state.reflecting_face = Geo_State_t::NONE;
    state.face            = Geo_State_t::NONE;
    state.region          = Geo_State_t::NONE;
    std::fill(state.exiting_level.begin(), state.exiting_level.end(), 0);

    // find the object that this point is in
    int id = d_layout[find_object(r, state)];
    CHECK(d_objects[id]);

    // Transformed coordinates.
    Space_Vector tr = transform(r, state);

    // initialize the state by diving into the object
    d_objects[id]->initialize(tr, state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Track distance to boundary.
 */
template<class T>
void RTK_Array<T>::distance_to_boundary(const Space_Vector &r,
                                        const Space_Vector &omega,
                                        Geo_State_t        &state) const
{
    REQUIRE(d_completed);

    // clear reflecting face indicator
    state.reflecting_face = Geo_State_t::NONE;

    // Transformed coordinates.
    Space_Vector tr = transform(r, state);

    // dive through objects until we hit the pin-cell
    object(state)->distance_to_boundary(tr, omega, state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build volumes.
 *
 * \note The vessel is not included in the volume calculations (the vessel
 * does not have a cellid).
 */
template<class T>
typename RTK_Array<T>::Vec_Dbl
RTK_Array<T>::get_volumes() const
{
    Vec_Dbl volumes(d_total_cells, 0.0);

    // build the volumes
    build_volumes(volumes, 0);

#ifdef REMEMBER_ON
    using def::X; using def::Y; using def::Z;
    double total = std::accumulate(volumes.begin(), volumes.end(), 0.0);
    double ref   = pitch(X) * pitch(Y) * height();
    VALIDATE(profugus::soft_equiv(total, ref, 1.0e-14 * volumes.size()),
             "Added cell volumes = " << total << " are not equal to the "
             << "total volume of " << ref);
#endif
    return volumes;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the state of particle while in mid-track (at collision sites).
 */
template<class T>
void RTK_Array<T>::update_state(Geo_State_t &state) const
{
    REQUIRE(d_completed);

    // clear reflecting face indicator
    state.reflecting_face = Geo_State_t::NONE;

    // update the state of the current object if it has not escaped
    if (state.escaping_face == Geo_State_t::NONE)
        object(state)->update_state(state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a surface.
 */
template<class T>
void RTK_Array<T>::cross_surface(const Space_Vector &r,
                                 Geo_State_t        &state)
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(d_completed);

    // initialize exiting level flag
    std::fill(state.exiting_level.begin(), state.exiting_level.end(), 0);

    // determine the surface crossings in each level
    determine_boundary_crossings(state);

    // simply return if we are at an internal pin face as everything should be
    // correctly set
    if (state.exiting_face == Geo_State_t::INTERNAL)
    {
        CHECK(!state.exiting_level[d_level]);
        return;
    }

    // update all of the downstream levels state assuming we haven't left the
    // geometry
    else if (!state.exiting_level[d_level])
    {
        update_coordinates(r, state);
    }

    // otherwise, the particle has escaped the geometry, record the state
    else
    {
        CHECK(state.exiting_level[d_level]);
        CHECK(state.exiting_face > Geo_State_t::INTERNAL);

        int refl_face_index = state.exiting_face - Geo_State_t::MINUS_X;
        CHECK(refl_face_index >= 0 && refl_face_index < 6);

        // if this is a reflecting face, then reflect the particle and add 1
        // to all of the level coordinates
        if (d_reflect[refl_face_index])
        {
            // store the reflecting face
            state.reflecting_face = state.exiting_face;

            // set the face to none representing an external pin face
            state.face = Geo_State_t::NONE;

            // get correction for logical coordinates of objects in each level
            Logical_Array reflect(0, 0, 0);

            switch (state.reflecting_face)
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
                    VALIDATE(false, "Lost a reflecting face.");
            }

            // loop through levels and update the coordinates to bring them
            // back into the geometry
            for (int l = 0; l <= d_level; ++l)
            {
                state.level_coord[l] += reflect;
            }
        }

        // otherwise, the particle has left the geometry
        else
        {
            state.escaping_face = state.exiting_face;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the object a point is in.
 */
template<class T>
int RTK_Array<T>::find_object(const Space_Vector &r,
                              Geo_State_t        &state) const
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(r[X] >= d_x.front()); REQUIRE(r[X] <= d_x.back());
    REQUIRE(r[Y] >= d_y.front()); REQUIRE(r[Y] <= d_y.back());
    REQUIRE(r[Z] >= d_z.front()); REQUIRE(r[Z] <= d_z.back());

    // iterators
    int                     i,j,k;
    Vec_Dbl::const_iterator itr, jtr, ktr;

    // find the logical indices of the object in the array
    itr = std::lower_bound(d_x.begin(), d_x.end(), r[X]);
    i   = itr - d_x.begin() - 1;
    jtr = std::lower_bound(d_y.begin(), d_y.end(), r[Y]);
    j   = jtr - d_y.begin() - 1;
    ktr = std::lower_bound(d_z.begin(), d_z.end(), r[Z]);
    k   = ktr - d_z.begin() - 1;

    // check for particles on the low face of the array
    if (r[X] == d_x[0])
    {
        CHECK(itr == d_x.begin());

        // reset index
        i = 0;
    }
    if (r[Y] == d_y[0])
    {
        CHECK(jtr == d_y.begin());

        // reset index
        j = 0;
    }
    if (r[Z] == d_z[0])
    {
        CHECK(ktr == d_z.begin());

        // reset index
        k = 0;
    }

    CHECK(i >= 0 && i < d_N[X]);
    CHECK(j >= 0 && j < d_N[Y]);
    CHECK(k >= 0 && k < d_N[Z]);

    ENSURE(d_x[i] <= r[X] && d_x[i+1] >= r[X]);
    ENSURE(d_y[j] <= r[Y] && d_y[j+1] >= r[Y]);
    ENSURE(d_z[k] <= r[Z] && d_z[k+1] >= r[Z]);
    state.level_coord[d_level][X] = i;
    state.level_coord[d_level][Y] = j;
    state.level_coord[d_level][Z] = k;

    return index(i, j, k);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the object on a boundary face.
 */
template<class T>
int RTK_Array<T>::find_object_on_boundary(const Space_Vector &r,
                                          int                 face,
                                          int                 face_type,
                                          Geo_State_t        &state) const
{
    using def::X; using def::Y; using def::Z;

    // iterators
    int                     i = 0, j = 0, k = 0;
    Vec_Dbl::const_iterator itr, jtr, ktr;

    // find the logical indices of the object in the array; do not search the
    // face dimension as this is known
    if (face_type != X)
    {
        itr = std::lower_bound(d_x.begin(), d_x.end(), r[X]);
        i   = itr - d_x.begin() - 1;
        CHECK(d_x[i] <= r[X] && d_x[i+1] >= r[X]);
    }
    else
    {
        if (face == Geo_State_t::PLUS_X)
            i = d_N[X] - 1;
    }

    if (face_type != Y)
    {
        jtr = std::lower_bound(d_y.begin(), d_y.end(), r[Y]);
        j   = jtr - d_y.begin() - 1;
        CHECK(d_y[j] <= r[Y] && d_y[j+1] >= r[Y]);
    }
    else
    {
        if (face == Geo_State_t::PLUS_Y)
            j = d_N[Y] - 1;
    }

    if (face_type != Z)
    {
        ktr = std::lower_bound(d_z.begin(), d_z.end(), r[Z]);
        k   = ktr - d_z.begin() - 1;
        CHECK(d_z[k] <= r[Z] && d_z[k+1] >= r[Z]);
    }
    else
    {
        if (face == Geo_State_t::PLUS_Z)
            k = d_N[Z] - 1;
    }

    ENSURE(i >= 0 && i < d_N[X]);
    ENSURE(j >= 0 && j < d_N[Y]);
    ENSURE(k >= 0 && k < d_N[Z]);

    state.level_coord[d_level][X] = i;
    state.level_coord[d_level][Y] = j;
    state.level_coord[d_level][Z] = k;

    return index(i, j, k);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output the array heirarchy for diagnostics.
 */
template<class T>
void RTK_Array<T>::output(std::ostream &out) const
{
    using std::endl; using std::setw; using std::setfill; using std::left;
    using std::right; using std::fixed; using std::setprecision;
    using def::I; using def::J; using def::K;

    // write basic header info
    out << "****************************************"
        << "****************************************" << endl;
    out << endl;
    out << left;

    out << "Number of levels   = " << setw(3) << d_level + 1 << endl;
    out << "Array map at level = " << setw(3) << d_level << endl;
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

    for (int n = 0; n < d_objects.size(); ++n)
    {
        if (d_objects[n])
        {
            d_objects[n]->output(out, d_level, n);
        }
    }
    out << endl;
    out << right;
    out << "****************************************"
        << "****************************************" << endl;
}

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the array level.
 */
template<class T>
int RTK_Array<T>::calc_level()
{
    return 1 + T::calc_level();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Total number of regions
 */
template<class T>
int RTK_Array<T>::num_regions() const
{
    int regions = 0;
    for (int i=0; i < size(def::X); ++i)
    {
        for (int j=0; j < size(def::Y); ++j)
        {
            for (int k = 0; k < size(def::Z); ++k)
            {
                regions += object(i,j,k).num_regions();
            }
        }
    }

    return regions;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Determine boundary crossings at each level starting at the lowest.
 */
template<class T>
void RTK_Array<T>::determine_boundary_crossings(Geo_State_t &state)
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(d_level > 0);

    // dive into the object and see if it crosses a boundary
    object(state)->determine_boundary_crossings(state);

    // process particles that cross an array boundary on the previous level,
    // they may cross a boundary at this level as well
    if (state.exiting_level[d_level - 1])
    {
        switch (state.exiting_face)
        {
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
                VALIDATE(false,
                    "Not at a valid array surface crossing.");
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the array coordinates in all levels where the particle has
 * crossed a face.
 */
template<class T>
void RTK_Array<T>::update_coordinates(const Space_Vector &r,
                                      Geo_State_t        &state)
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(d_level > 0);

    // transform coordinates into object's coordinate system
    Space_Vector tr = transform(r, state);

    // if the child object has exited the level, then update the coordinates
    if (state.exiting_level[d_level - 1])
    {
        switch (state.exiting_face)
        {
            // update the coordinates of the object across the given face
            case Geo_State_t::MINUS_X:
                object(state)->find_object_on_boundary(
                    tr, Geo_State_t::PLUS_X, X, state);
                break;
            case Geo_State_t::PLUS_X:
                object(state)->find_object_on_boundary(
                    tr, Geo_State_t::MINUS_X, X, state);
                break;
            case Geo_State_t::MINUS_Y:
                object(state)->find_object_on_boundary(
                    tr, Geo_State_t::PLUS_Y, Y, state);
                break;
            case Geo_State_t::PLUS_Y:
                object(state)->find_object_on_boundary(
                    tr, Geo_State_t::MINUS_Y, Y, state);
                break;
            case Geo_State_t::MINUS_Z:
                object(state)->find_object_on_boundary(
                    tr, Geo_State_t::PLUS_Z, Z, state);
                break;
            case Geo_State_t::PLUS_Z:
                object(state)->find_object_on_boundary(
                    tr, Geo_State_t::MINUS_Z, Z, state);
                break;
        }
    }

    // go to child object and recursively update all coordinates
    object(state)->update_coordinates(tr, state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Recursively count the number of cells in each object.
 */
template<class T>
int RTK_Array<T>::count_cells()
{
    using def::I; using def::J; using def::K;

    // initialize total cell count
    d_total_cells = 0;

    // array index
    int n = 0;

    // loop through objects and count the number of cells in each one
    for (int k = 0; k < d_N[K]; ++k)
    {
        for (int j = 0; j < d_N[J]; ++j)
        {
            for (int i = 0; i < d_N[I]; ++i)
            {
                // calculate the array index
                n = index(i, j, k);
                CHECK(n < d_num_cells.size());

                // count the cells in the object
                d_num_cells[n] = cell_count_dispatch(i, j, k);

                // calculate the total running sum
                d_total_cells += d_num_cells[n];
            }
        }
    }

    ENSURE(d_total_cells > 0);
    return d_total_cells;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add vessel to objects in the array.
 *
 * \param xoff distance from the origin of the vessel to the \b left side of
 * the object
 *
 * \param yoff distance from the origin of the vessel to the \b bottom of the
 * of the object
 */
template<class T>
void RTK_Array<T>::add_vessel(double R0,
                              double R1,
                              double xoff,
                              double yoff,
                              int    id)
{
    using def::X; using def::Y; using def::Z;

    // near and far points for the current object
    double near[2] = {0.0}, far[2] = {0.0};

    // near and far radii for the current object
    double nearR2 = 0.0, farR2 = 0.0;

    // vessel radii squared
    double R02 = R0*R0;
    double R12 = R1*R1;

    // running x,y offsets
    double rox = xoff, roy = yoff;

    // loop through objects in this array
    for (int j = 0; j < d_N[Y]; ++j)
    {
        // reset running x-offset
        rox = xoff;

        for (int i = 0; i < d_N[X]; ++i)
        {
            // check to see if the vessel bisects this array object
            near[X] = rox;
            far[X]  = rox + dx(i);
            near[Y] = roy;
            far[Y]  = roy + dy(j);

            if (xoff < 0.0)
            {
                near[X] = far[X];
                far[X]  = rox;
            }

            if (yoff < 0.0)
            {
                near[Y] = far[Y];
                far[Y]  = roy;
            }

            // the "=" on the check is for -x to +x on centerlines where x ==
            // x
            CHECK(std::fabs(near[X]) < std::fabs(far[X]) ||
                   soft_equiv(std::fabs(near[X]), std::fabs(far[X])));
            CHECK(std::fabs(near[Y]) < std::fabs(far[Y]) ||
                   soft_equiv(std::fabs(near[Y]), std::fabs(far[Y])));

            // near and far radii bounding the cell object
            nearR2 = near[X]*near[X] + near[Y]*near[Y];
            farR2  = far[X]*far[X] + far[Y]*far[Y];

            // now check to see if one of the vessel radii bisects the
            // cell
            if (!(R02 > farR2 || R12 < nearR2))
            {
                if (R02 >= nearR2 || R12 < farR2)
                {
                    // dive into this object and add the vessel
                    for (int k = 0; k < d_N[Z]; ++k)
                    {
                        add_vessel_to_object(i, j, k, R0, R1, rox, roy, id);
                    }
                }
            }

            // update the running xoffset
            rox += dx(i);

            CHECK(rox - xoff < pitch(X) || soft_equiv(rox - xoff, pitch(X)));
        }

        // update the running y-offset
        roy += dy(j);
    }

    ENSURE(soft_equiv(rox - xoff, pitch(X)));
    ENSURE(soft_equiv(roy - yoff, pitch(Y)));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build volumes.
 */
template<class T>
void RTK_Array<T>::build_volumes(Vec_Dbl &v,
                                 int      offset) const
{
    using def::X; using def::Y; using def::Z;

    // recursively dive into objects and calculate the volumes
    for (int k = 0; k < d_N[Z]; ++k)
    {
        for (int j = 0; j < d_N[Y]; ++j)
        {
            for (int i = 0; i < d_N[X]; ++i)
            {
                int off = d_Nc_offset[index(i, j, k)] + offset;
                object(i,j,k).build_volumes(v, off);
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output the array heirarchy for diagnostics.
 */
template<class T>
void RTK_Array<T>::output(std::ostream &out,
                          int           level,
                          int           obj_id) const
{
    using std::endl; using std::setw; using std::setfill; using std::left;
    using std::right; using std::fixed; using std::setprecision;
    using def::I; using def::J; using def::K;

    out << "----------------------------------------"
        << "----------------------------------------" << endl;
    out << endl;
    out << left;

    out << "Output array " << setw(3) << obj_id << " from level "
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
            d_objects[n]->output(out, d_level, n);
        }
    }
    out << endl;
    out << right;
    out << "----------------------------------------"
        << "----------------------------------------" << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return bounding box for given cell
 */
template<class T>
void RTK_Array<T>::get_cell_extents(geometry::cell_type cell,
                                    Space_Vector       &lower,
                                    Space_Vector       &upper) const
{
    REQUIRE( cell < d_total_cells );

    using def::I;
    using def::J;
    using def::K;

    // Locate cell in array
    int elem = std::upper_bound(d_Nc_offset.begin(),d_Nc_offset.end(),cell)
               - d_Nc_offset.begin() - 1;
    ENSURE( elem < d_Nc_offset.size() );

    // Get array element at this location and compute it's bounding box
    CHECK( elem < d_layout.size() );
    CHECK( d_layout[elem] < d_objects.size() );
    auto obj = d_objects[d_layout[elem]];
    ENSURE( obj );

    // Get local bounding box for this element
    CHECK( d_Nc_offset[elem] <= cell );
    obj->get_cell_extents(cell - d_Nc_offset[elem], lower, upper);

    // Offset by lower corner of object in its local frame
    auto obj_corner = obj->corner();
    for( auto ijk : {I, J, K} )
    {
        lower[ijk] -= obj_corner[ijk];
        upper[ijk] -= obj_corner[ijk];
    }

    // Compute ijk index of element
    int k = elem / (d_N[I] * d_N[J]);
    elem = elem % (d_N[I] * d_N[J]);
    int j = elem / d_N[I];
    int i = elem % d_N[I];

    // Offset by lower corner of object in this frame
    Space_Vector corner = {d_x[i], d_y[j], d_z[k]};
    lower += corner;
    upper += corner;
}

} // end namespace profugus

#endif // geometry_RTK_Array_t_hh

//---------------------------------------------------------------------------//
//                   end of geometry/RTK_Array.t.hh
//---------------------------------------------------------------------------//
