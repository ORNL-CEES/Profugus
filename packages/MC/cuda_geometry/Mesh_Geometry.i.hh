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

    ENSURE(state.ijk.i >= -1 && state.ijk.i <= d_mesh.num_cells_along(I));
    ENSURE(state.ijk.j >= -1 && state.ijk.j <= d_mesh.num_cells_along(J));
    ENSURE(state.ijk.k >= -1 && state.ijk.i <= d_mesh.num_cells_along(K));

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

    REQUIRE(soft_equiv(vector_magnitude(state.d_dir), 1.0, 1.e-5));

    const double * edges_x(d_mesh.edges(I));
    const double * edges_y(d_mesh.edges(J));
    const double * edges_z(d_mesh.edges(K));
    cuda::Coordinates extents = {d_mesh.num_cells_along(I),
                                 d_mesh.num_cells_along(J),
                                 d_mesh.num_cells_along(K)};

    // min distance to boundary; initializing test_dist to a large number before
    // each surface check implicitly handles the case where omega[dir] == 0.0
    double test_dist = cuda::constants::max_double;
    int    test_next = state.ijk.i;

    // initialize the next surface
    state.next_ijk = state.ijk;

    // unrolled check

    // X SURFACE
    if (state.d_dir.x > 0.0 && state.ijk.i < extents.i)
    {
        test_dist = (edges_x[state.ijk.i+1] - state.d_r.x) / state.d_dir.x;
        test_next = state.ijk.i + 1;
    }
    else if (state.d_dir.x < 0.0 && state.ijk.i > -1)
    {
        test_dist = (edges_x[state.ijk.i] - state.d_r.x) / state.d_dir.x;
        test_next = state.ijk.i - 1;
    }

    // initialize to x distance
    state.next_dist   = test_dist;
    state.next_ijk.i = test_next;

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = cuda::constants::max_double;

    // Y SURFACE
    if (state.d_dir.y > 0.0 && state.ijk.j < extents.j)
    {
        test_dist = (edges_y[state.ijk.j+1] - state.d_r.y) / state.d_dir.y;
        test_next = state.ijk.j + 1;
    }
    else if (state.d_dir.y < 0.0 && state.ijk.j > -1)
    {
        test_dist = (edges_y[state.ijk.j] - state.d_r.y) / state.d_dir.y;
        test_next = state.ijk.j - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist   = test_dist;
        state.next_ijk.i = state.ijk.i;
        state.next_ijk.j = test_next;
    }

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = cuda::constants::max_double;

    // Z SURFACE
    if (state.d_dir.z > 0.0 && state.ijk.k < extents.k)
    {
        test_dist = (edges_z[state.ijk.k+1] - state.d_r.z) / state.d_dir.z;
        test_next = state.ijk.k + 1;
    }
    else if (state.d_dir.z < 0.0 && state.ijk.k > -1)
    {
        test_dist = (edges_z[state.ijk.k] - state.d_r.z) / state.d_dir.z;
        test_next = state.ijk.k - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist   = test_dist;
        state.next_ijk.i = state.ijk.i;
        state.next_ijk.j = state.ijk.j;
        state.next_ijk.k = test_next;
    }

    ENSURE(state.next_dist >= 0.);
    return state.next_dist;
}

} // end namespace cuda_profugus

#endif // MC_cuda_geometry_Mesh_Geometry_i_hh

//---------------------------------------------------------------------------//
//                 end of Mesh_Geometry.i.hh
//---------------------------------------------------------------------------//
