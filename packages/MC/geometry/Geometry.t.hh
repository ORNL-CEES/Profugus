//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Geometry.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 25 10:02:33 2011
 * \brief  RTK_Geometry template member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_RTK_Geometry_t_hh
#define geometry_RTK_Geometry_t_hh

#include "RTK_Geometry.hh"
#include "RTK_Functions.hh"

namespace denovo
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Array>
RTK_Geometry<Array>::RTK_Geometry(SP_Array array)
    : d_array(array)
    , d_volumes(d_array->num_cells(), 0.0)
    , d_level(d_array->level())
    , d_modular( false )
{
    d_array->get_extents(d_lower, d_upper);
    Ensure (d_array);
    Ensure (d_lower <= d_upper);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Array>
RTK_Geometry<Array>::RTK_Geometry(SP_Array array, bool symmetric)
    : d_array(array)
    , d_volumes(d_array->num_cells(), 0.0)
    , d_mapped_cells(d_array->num_cells(), -1)
    , d_level(d_array->level())
    , d_modular( false )
{
    d_array->get_extents(d_lower, d_upper);
    Ensure (d_array);
    Ensure (d_lower <= d_upper);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*! \brief Initialize a track.
 *
 * This also checks for whether the starting position is outside the geometry
 * (possible for sources and especially meshing).
 */
template<class Array>
void RTK_Geometry<Array>::initialize(const Space_Vector &r,
                                     const Space_Vector &direction,
                                     Geo_State_t        &state)
{
    Require (d_array);

    // add position and direction to the state
    state.d_r   = r;
    state.d_dir = direction;

    // normalize the direction
    vector_normalize(state.d_dir);

    // Check whether the point is in the geometry. We have to preserve the
    // escaping face because the underlying elements reset it by default.
    using def::X; using def::Y; using def::Z;
    int escaping = Geo_State_t::NONE;
    if (   (r[X] < d_lower[X]) || (r[X] > d_upper[X])
        || (r[Y] < d_lower[Y]) || (r[Y] > d_upper[Y])
        || (r[Z] < d_lower[Z]) || (r[Z] > d_upper[Z]))
    {
        const Space_Vector& omega = state.d_dir;
        Validate(
               (omega[X] < 0 ? r[X] > d_lower[X] : r[X] < d_upper[X])
            || (omega[Y] < 0 ? r[Y] > d_lower[Y] : r[Y] < d_upper[Y])
            || (omega[Z] < 0 ? r[Z] > d_lower[Z] : r[Z] < d_upper[Z]),
            "Particle started outside the geometry at " << r << " in "
            "direction " << omega);

        state.escaping_face = Geo_State_t::NONE;
        move_from_outside(d_lower, d_upper, state);
        escaping = state.escaping_face;
    }

    // initialize the array with the current position
    d_array->initialize(state.d_r, state);

    // We have to reset the escaping parameter because "initialize" will
    // overwrite it.
    if (escaping != Geo_State_t::NONE)
    {
        state.escaping_face    = escaping;
        state.exiting_face     = escaping;
        state.exiting_level    = escaping; // assign to all levels
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reflect the direction at a reflecting surface with outgoing normal
 * \f$\hat{\mathbf{n}}\f$.
 *
 * The reflected angle off of a surface with outgoing normal
 * \f$\hat{\mathbf{n}}\f$ is
 * \f[
   \hat{\Omega}_{m'}=\hat{\Omega}_m - 2\hat{\mathbf{n}}(\hat{Omega}\cdot
   \hat{\mathbf{n}})\:.
 * \f]
 *
 * \return true if direction was reflected; false if not at a reflecting
 * surface
 */
template<class Array>
bool RTK_Geometry<Array>::reflect(Geo_State_t &state)
{
    using def::X; using def::Y; using def::Z;

    Require (nemesis::soft_equiv(
                 vector_magnitude(state.d_dir), 1.0, 1.0e-6));

    // get the outward normal
    Space_Vector n = normal(state);

    // calculate the dot-product of the incoming angle and outward normal
    double dot = state.d_dir[X]*n[X] + state.d_dir[Y]*n[Y] + state.d_dir[Z]*n[Z];

    // if the dot-product != 0 then calculate the reflected angle
    if (dot != 0.0)
    {
        state.d_dir[X] -= 2.0 * n[X] * dot;
        state.d_dir[Y] -= 2.0 * n[Y] * dot;
        state.d_dir[Z] -= 2.0 * n[Z] * dot;

        Ensure (nemesis::soft_equiv(
                    vector_magnitude(state.d_dir), 1.0, 1.0e-6));
        return true;
    }

    // we were not at a reflecting surface
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the outward normal for points on a face.
 */
template<class Array>
typename RTK_Geometry<Array>::Space_Vector
RTK_Geometry<Array>::normal(const Geo_State_t &state) const
{
    // query is based on position of particle on a face; otherwise we return a
    // zero vector
    switch (state.exiting_face)
    {
        case Geo_State_t::MINUS_X:
            return Space_Vector(-1.0, 0.0, 0.0);

        case Geo_State_t::PLUS_X:
            return Space_Vector(1.0, 0.0, 0.0);

        case Geo_State_t::MINUS_Y:
            return Space_Vector(0.0, -1.0, 0.0);

        case Geo_State_t::PLUS_Y:
            return Space_Vector(0.0, 1.0, 0.0);

        case Geo_State_t::MINUS_Z:
            return Space_Vector(0.0, 0.0, -1.0);

        case Geo_State_t::PLUS_Z:
            return Space_Vector(0.0, 0.0, 1.0);

        default:
            break;
    }

    // return 0 vector if not on an exiting face
    return Space_Vector(0.0, 0.0, 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sets mapped cells
 */
template<class Array>
void RTK_Geometry<Array>::set_mapped_cells()
{
    d_array->set_mapped_cells(
        0, 0, d_mapped_cells.begin(), d_mapped_cells.begin());

} // end set_mapped_cells

} // end namespace denovo

#endif // geometry_RTK_Geometry_t_hh

//---------------------------------------------------------------------------//
//                   end of geometry/RTK_Geometry.t.hh
//---------------------------------------------------------------------------//
