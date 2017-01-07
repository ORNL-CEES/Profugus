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
void RTK_Array<T>::initialize(const Space_Vector &r,
                              Geo_State_t        &state) const
{
    // initialize state
    state.escaping_face   = Geo_State_t::NONE;
    state.reflecting_face = Geo_State_t::NONE;
    state.face            = Geo_State_t::NONE;
    state.region          = Geo_State_t::NONE;
    state.exiting_level   = {0, 0, 0};

    // find the object that this point is in
    int id = d_layout[find_object(r, state)];

    // Transformed coordinates.
    Space_Vector tr = transform(r, state);

    // initialize the state by diving into the object
    d_objects[id].initialize(tr, state);
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
    Geo_State_t        &state) const
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
    state.level_coord[d_level][X] = i;
    state.level_coord[d_level][Y] = j;
    state.level_coord[d_level][Z] = k;

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
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Transform into object coordinate system.
 */
template<class T>
__device__
typename RTK_Array<T>::Space_Vector
RTK_Array<T>::transform(const Space_Vector &r,
                        const Geo_State_t  &state) const
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
__device__
const typename RTK_Array<T>::Object_t&
RTK_Array<T>::object(const Geo_State_t &state) const
{
    using def::X; using def::Y; using def::Z;
    return d_objects[d_layout[index(state.level_coord[d_level][X],
                                    state.level_coord[d_level][Y],
                                    state.level_coord[d_level][Z])]];
}

//---------------------------------------------------------------------------//
// SPECIALIZATIONS ON RTK_CELL
//---------------------------------------------------------------------------//

template<>
__device__
inline RTK_Array<RTK_Cell>::Space_Vector
RTK_Array<RTK_Cell>::transform(const Space_Vector &r,
                               const Geo_State_t  &state) const
{
    using def::X; using def::Y; using def::Z;
    using cuda::utility::soft_equiv;
    using std::fabs;

    DEVICE_REQUIRE(d_level == 0);

    // Transformed coordinates and lower/upper extents
    Space_Vector tr, lower, upper;

    // get the extents (we only need the lower)
    object(state).get_extents(lower, upper);
    DEVICE_CHECK(lower[X] < 0.0);
    DEVICE_CHECK(lower[Y] < 0.0);

    // transform the coordinates to the pin cell
    tr[X] = (r[X] - d_x[state.level_coord[0][X]]) + lower[X];
    tr[Y] = (r[Y] - d_y[state.level_coord[0][Y]]) + lower[Y];
    tr[Z] = (r[Z] - d_z[state.level_coord[0][Z]]);

    DEVICE_ENSURE((tr[X] > lower[X] && tr[X] < upper[X]) ||
                  (soft_equiv(tr[X], lower[X], 1.0e-6 * fabs(lower[X])) ||
                   soft_equiv(tr[X], upper[X], 1.0e-6 * upper[X])));
    DEVICE_ENSURE((tr[Y] > lower[Y] && tr[Y] < upper[Y]) ||
                  (soft_equiv(tr[Y], lower[Y], 1.0e-6 * fabs(lower[Y])) ||
                   soft_equiv(tr[Y], upper[Y], 1.0e-6 * upper[Y])));
    return tr;
}

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_Array_i_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.i.cuh
//---------------------------------------------------------------------------//
