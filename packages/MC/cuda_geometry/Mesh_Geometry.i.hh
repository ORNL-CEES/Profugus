//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_geometry/Mesh_Geometry.cc
 * \author Thomas M. Evans and Seth Johnson
 * \date   Mon Jul 21 17:56:40 2014
 * \brief  Mesh_Geometry member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_Mesh_Geometry_i_hh
#define MC_cuda_geometry_Mesh_Geometry_i_hh

#include "Mesh_Geometry.hh"

#include "utils/Definitions.hh"
#include "cuda_utils/Constants.hh"
#include "cuda_utils/CudaDBC.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize particle
 *
 * This finds the closest ijk coords on each axis to the location. A particle
 * can be born "outside" and have ijk extents that are outside [0,N) .
 */
__device__ void Mesh_Geometry::initialize(
        const Space_Vector& r,
        const Space_Vector& direction,
        Geo_State_t       & state) const
{
    using cuda::utility::soft_equiv;
    using cuda::utility::vector_normalize;
    using def::I; using def::J; using def::K;

    // Set struct attributes
    state.d_r = r;
    state.d_dir = direction;

    vector_normalize(state.d_dir);

    update_state(state);

    DEVICE_ENSURE(state.ijk[I] >= -1 &&
                  state.ijk[I] <= d_mesh.num_cells_along(I));
    DEVICE_ENSURE(state.ijk[J] >= -1 &&
                  state.ijk[J] <= d_mesh.num_cells_along(J));
    DEVICE_ENSURE(state.ijk[K] >= -1 &&
                  state.ijk[K] <= d_mesh.num_cells_along(K));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate distance to the next cell
 */
__device__
double Mesh_Geometry::distance_to_boundary(Geo_State_t& state) const
{
    using def::I; using def::J; using def::K;
    using cuda::utility::soft_equiv;
    using cuda::utility::vector_magnitude;

    DEVICE_REQUIRE(soft_equiv(vector_magnitude(state.d_dir), 1.0, 1.e-5));

    const Double_View& edges_x = d_mesh.edges(I);
    const Double_View& edges_y = d_mesh.edges(J);
    const Double_View& edges_z = d_mesh.edges(K);

    cuda_utils::Coordinates extents = {d_mesh.num_cells_along(I),
                                       d_mesh.num_cells_along(J),
                                       d_mesh.num_cells_along(K)};

    // min distance to boundary; initializing test_dist to a large number before
    // each surface check implicitly handles the case where omega[dir] == 0.0
    double test_dist = cuda::constants::max_double;
    int    test_next = state.ijk[I];

    // initialize the next surface
    state.next_ijk = state.ijk;

    // unrolled check

    // X SURFACE
    if (state.d_dir[I] > 0.0 && state.ijk[I] < extents[I])
    {
        test_dist = (edges_x[state.ijk[I]+1] - state.d_r[I]) / state.d_dir[I];
        test_next = state.ijk[I] + 1;
    }
    else if (state.d_dir[I] < 0.0 && state.ijk[I] > -1)
    {
        test_dist = (edges_x[state.ijk[I]] - state.d_r[I]) / state.d_dir[I];
        test_next = state.ijk[I] - 1;
    }

    // initialize to x distance
    state.next_dist   = test_dist;
    state.next_ijk[I] = test_next;

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = cuda::constants::max_double;

    // Y SURFACE
    if (state.d_dir[J] > 0.0 && state.ijk[J] < extents[J])
    {
        test_dist = (edges_y[state.ijk[J]+1] - state.d_r[J]) / state.d_dir[J];
        test_next = state.ijk[J] + 1;
    }
    else if (state.d_dir[J] < 0.0 && state.ijk[J] > -1)
    {
        test_dist = (edges_y[state.ijk[J]] - state.d_r[J]) / state.d_dir[J];
        test_next = state.ijk[J] - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist   = test_dist;
        state.next_ijk[I] = state.ijk[I];
        state.next_ijk[J] = test_next;
    }

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = cuda::constants::max_double;

    // Z SURFACE
    if (state.d_dir[K] > 0.0 && state.ijk[K] < extents[K])
    {
        test_dist = (edges_z[state.ijk[K]+1] - state.d_r[K]) / state.d_dir[K];
        test_next = state.ijk[K] + 1;
    }
    else if (state.d_dir[K] < 0.0 && state.ijk[K] > -1)
    {
        test_dist = (edges_z[state.ijk[K]] - state.d_r[K]) / state.d_dir[K];
        test_next = state.ijk[K] - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist  = test_dist;
        state.next_ijk[I] = state.ijk[I];
        state.next_ijk[J] = state.ijk[J];
        state.next_ijk[K] = test_next;
    }

    DEVICE_ENSURE(state.next_dist >= 0.);
    return state.next_dist;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reflect the direction on a reflecting surface
 */
__device__
bool Mesh_Geometry::reflect(Geo_State_t& state) const
{
    using def::I; using def::J; using def::K;

    using cuda::utility::soft_equiv;
    using cuda::utility::vector_magnitude;
    DEVICE_REQUIRE(soft_equiv(vector_magnitude(state.d_dir), 1.0, 1.0e-6));

    // If we're not in a reflecting state, return
    if( state.reflecting_face == Geo_State_t::NONE )
        return false;

    // get the outward normal
    Space_Vector n = normal(state);

    // calculate the dot-product of the incoming angle and outward normal
    double dot = state.d_dir[I]*n[I] +
                 state.d_dir[J]*n[J] +
                 state.d_dir[K]*n[K];

    DEVICE_CHECK( dot != 0.0 );

    state.d_dir[I] -= 2.0 * n[I] * dot;
    state.d_dir[J] -= 2.0 * n[J] * dot;
    state.d_dir[K] -= 2.0 * n[K] * dot;

    DEVICE_ENSURE(soft_equiv(vector_magnitude(state.d_dir), 1.0, 1.0e-6));

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the outward normal at the location dictated by the state.
 */
__device__ cuda_utils::Space_Vector
Mesh_Geometry::normal(const Geo_State_t& state) const
{
    // Choose normal based on exiting face
    if( state.exiting_face == Geo_State_t::MINUS_X )
    {
        return {-1.0, 0.0, 0.0};
    }
    else if( state.exiting_face == Geo_State_t::PLUS_X )
    {
        return {1.0, 0.0, 0.0};
    }
    else if( state.exiting_face == Geo_State_t::MINUS_Y )
    {
        return {0.0, -1.0, 0.0};
    }
    else if( state.exiting_face == Geo_State_t::PLUS_Y )
    {
        return {0.0, 1.0, 0.0};
    }
    else if( state.exiting_face == Geo_State_t::MINUS_Z )
    {
        return {0.0, 0.0, -1.0};
    }
    else if( state.exiting_face == Geo_State_t::PLUS_Z )
    {
        return {0.0, 0.0, 1.0};
    }

    // We weren't on a boundary
    return {0.0, 0.0, 0.0};
}

} // end namespace cuda_profugus

#endif // MC_cuda_geometry_Mesh_Geometry_i_hh

//---------------------------------------------------------------------------//
//                 end of Mesh_Geometry.i.hh
//---------------------------------------------------------------------------//
