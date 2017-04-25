//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Geometry.i.cuh
 * \author Tom Evans
 * \date   Mon Jan 30 00:14:27 2017
 * \brief  RTK_Geometry CUDA device class definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Geometry_i_cuh
#define MC_cuda_rtk_RTK_Geometry_i_cuh

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a track.
 */
__device__
inline void RTK_Geometry::initialize(
    const Space_Vector &r,
    const Space_Vector &direction,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    // add position and direction to the state
    state_vector.pos(pid) = r;
    state_vector.dir(pid) = direction;

    // normalize the direction
    cuda::utility::vector_normalize(state_vector.dir(pid));

    // initialize the array with the current position
    d_array.initialize(state_vector.pos(pid), state_vector, pid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the state with respect to outer geometry boundary.
 */
__device__
inline profugus::geometry::Boundary_State
RTK_Geometry::boundary_state(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    if (state_vector.escaping_face(pid) != Geo_State_t::NONE)
    {
        // if the particle has escaped indicate that the particle
        // is outside the geometry
        DEVICE_CHECK(state_vector.exiting_level(pid,d_level) ||
                     state_vector.next_face(pid) == Geo_State_t::NONE);
        return profugus::geometry::OUTSIDE;
    }
    else if (state_vector.reflecting_face(pid) != Geo_State_t::NONE)
    {
        // test of reflection on the given face
        return profugus::geometry::REFLECT;
    }
    return profugus::geometry::INSIDE;
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
__device__
inline bool RTK_Geometry::reflect(
    Geo_State_Vector_t &state_vector, int pid) const
{
    using def::X; using def::Y; using def::Z;

    DEVICE_REQUIRE(
        cuda::utility::soft_equiv(
            cuda::utility::vector_magnitude(state_vector.dir(pid)), 1.0, 1.0e-6));

    // get the outward normal
    Space_Vector n = normal(state_vector, pid);

    // calculate the dot-product of the incoming angle and outward normal
    double dot = state_vector.dir(pid)[X]*n[X] + state_vector.dir(pid)[Y]*n[Y] +
                 state_vector.dir(pid)[Z]*n[Z];

    // if the dot-product != 0 then calculate the reflected angle
    if (dot != 0.0)
    {
        state_vector.dir(pid)[X] -= 2.0 * n[X] * dot;
        state_vector.dir(pid)[Y] -= 2.0 * n[Y] * dot;
        state_vector.dir(pid)[Z] -= 2.0 * n[Z] * dot;

        DEVICE_ENSURE(
            cuda::utility::soft_equiv(
                cuda::utility::vector_magnitude(state_vector.dir(pid)),
                1.0, 1.0e-6));
        return true;
    }

    // we were not at a reflecting surface
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the outward normal for points on a face.
 */
__device__
inline RTK_Geometry::Space_Vector
RTK_Geometry::normal(
    const Geo_State_Vector_t &state_vector, int pid) const
{
    // query is based on position of particle on a face; otherwise we return a
    // zero vector
    switch (state_vector.exiting_face(pid))
    {
        case Geo_State_t::MINUS_X:
            return Space_Vector{-1.0, 0.0, 0.0};

        case Geo_State_t::PLUS_X:
            return Space_Vector{1.0, 0.0, 0.0};

        case Geo_State_t::MINUS_Y:
            return Space_Vector{0.0, -1.0, 0.0};

        case Geo_State_t::PLUS_Y:
            return Space_Vector{0.0, 1.0, 0.0};

        case Geo_State_t::MINUS_Z:
            return Space_Vector{0.0, 0.0, -1.0};

        case Geo_State_t::PLUS_Z:
            return Space_Vector{0.0, 0.0, 1.0};

        default:
            break;
    }

    // return 0 vector if not on an exiting face
    return Space_Vector{0.0, 0.0, 0.0};
}

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_Geometry_i_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Geometry.i.cuh
//---------------------------------------------------------------------------//
