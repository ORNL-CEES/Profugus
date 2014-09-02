//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Mesh_Geometry.cc
 * \author Thomas M. Evans and Seth Johnson
 * \date   Mon Jul 21 17:56:40 2014
 * \brief  Mesh_Geometry member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Mesh_Geometry.hh"

#include <limits>
#include <algorithm>

#include "harness/DBC.hh"
#include "utils/Container_Functions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Construct from cell edges
 */
Mesh_Geometry::Mesh_Geometry(
        const Vec_Dbl& x_edges,
        const Vec_Dbl& y_edges,
        const Vec_Dbl& z_edges)
    : d_mesh(x_edges, y_edges, z_edges)
{
    INSIST(d_mesh.dimension() == 3, "Only 3-D meshes are supported.");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set a field of material IDs (optional)
 */
void Mesh_Geometry::set_matids(SP_Vec_Int matids)
{
    REQUIRE(matids);
    REQUIRE(matids->size() == num_cells());
    d_materials = matids;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize particle
 *
 * This finds the closest ijk coords on each axis to the location. A particle
 * can be born "outside" and have ijk extents that are outside [0,N) .
 */
void Mesh_Geometry::initialize(
        const Space_Vector& r,
        const Space_Vector& direction,
        Geo_State_t       & state)
{
    using profugus::soft_equiv;
    using def::I; using def::J; using def::K;

    REQUIRE(soft_equiv(vector_magnitude(direction), 1.0, 1.e-5));

    // Set struct attributes
    state.pos = r;
    state.dir = direction;

    update_state(state);

    ENSURE(state.ijk[I] >= -1 && state.ijk[I] <= d_mesh.num_cells_along(I));
    ENSURE(state.ijk[J] >= -1 && state.ijk[J] <= d_mesh.num_cells_along(J));
    ENSURE(state.ijk[K] >= -1 && state.ijk[K] <= d_mesh.num_cells_along(K));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate distance to the next cell
 */
double Mesh_Geometry::distance_to_boundary(Geo_State_t& state)
{
    using def::I; using def::J; using def::K;
    using profugus::soft_equiv;

    REQUIRE(soft_equiv(vector_magnitude(state.dir), 1.0, 1.e-5));

    const Vec_Dbl& edges_x(d_mesh.edges(I));
    const Vec_Dbl& edges_y(d_mesh.edges(J));
    const Vec_Dbl& edges_z(d_mesh.edges(K));
    const Cartesian_Mesh::Dim_Vector& extents = d_mesh.extents();

    // min distance to boundary; initializing test_dist to a large number before
    // each surface check implicitly handles the case where omega[dir] == 0.0
    double test_dist = std::numeric_limits<double>::max();
    int    test_next = state.ijk[I];

    // initialize the next surface
    state.next_ijk = state.ijk;

    // unrolled check

    // X SURFACE
    if (state.dir[I] > 0.0 && state.ijk[I] < extents[I])
    {
        test_dist = (edges_x[state.ijk[I]+1] - state.pos[I]) / state.dir[I];
        test_next = state.ijk[I] + 1;
    }
    else if (state.dir[I] < 0.0 && state.ijk[I] > -1)
    {
        test_dist = (edges_x[state.ijk[I]] - state.pos[I]) / state.dir[I];
        test_next = state.ijk[I] - 1;
    }

    // initialize to x distance
    state.next_dist   = test_dist;
    state.next_ijk[I] = test_next;

    // reset the local dist-to-boundary to a large value to handle dir=0.0
    test_dist = std::numeric_limits<double>::max();

    // Y SURFACE
    if (state.dir[J] > 0.0 && state.ijk[J] < extents[J])
    {
        test_dist = (edges_y[state.ijk[J]+1] - state.pos[J]) / state.dir[J];
        test_next = state.ijk[J] + 1;
    }
    else if (state.dir[J] < 0.0 && state.ijk[J] > -1)
    {
        test_dist = (edges_y[state.ijk[J]] - state.pos[J]) / state.dir[J];
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
    test_dist = std::numeric_limits<double>::max();

    // Z SURFACE
    if (state.dir[K] > 0.0 && state.ijk[K] < extents[K])
    {
        test_dist = (edges_z[state.ijk[K]+1] - state.pos[K]) / state.dir[K];
        test_next = state.ijk[K] + 1;
    }
    else if (state.dir[K] < 0.0 && state.ijk[K] > -1)
    {
        test_dist = (edges_z[state.ijk[K]] - state.pos[K]) / state.dir[K];
        test_next = state.ijk[K] - 1;
    }

    // update running value of distance to boundary
    if (test_dist < state.next_dist)
    {
        state.next_dist   = test_dist;
        state.next_ijk[I] = state.ijk[I];
        state.next_ijk[J] = state.ijk[J];
        state.next_ijk[K] = test_next;
    }

    ENSURE(state.next_dist >= 0.);
    return state.next_dist;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the distance for moving the particle to an interior bound
 *
 * This function is provided so that if the particle misses the geometry
 * completely, you don't have to step through a bunch of non-existent cells in
 * your tally mesh.
 *
 * It should be called after initialize() to move the particle to the correct
 * inner boundary location and provide the distance that it moved.
 *
 * We don't do this in the "initialize" call because (e.g. in the Tally Mesh) we
 * want to actually keep track of the distance traveled outside this geometry's
 * domains.
 *
 * After calling this function, if the return value is numeric_limits::max(),
 * the next_ijk data is invalid. So after calling this and "move_to_surface",
 * which transports the particle to infinity, the state is invalid.
 *
 *
 * If you choose not to use the distance_to_interior function, the standard
 * distance tracking algorithm will work; although you can't use the shortcut of
 * breaking out when geo.cell_index == -1, because multiple distance_to_surface
 * calls will return surfaces outside the cartesian mesh domain.
 */
double Mesh_Geometry::distance_to_interior(Geo_State_t& state)
{
    using def::I; using def::J; using def::K;
    using def::X; using def::Y; using def::Z;
    using profugus::soft_equiv;

    // Initialize next distance to -1 for cases where particle begins on edges.
    state.next_dist = -1.;
    Geo_State_t::Face next_face = Geo_State_t::END_FACES;

    const Cartesian_Mesh::Dim_Vector& extents = d_mesh.extents();

    // First, check if we start outside, heading the wrong direction
    if (   (state.ijk[I] == -1         && state.dir[X] < 0.0)
           || (state.ijk[J] == -1         && state.dir[Y] < 0.0)
           || (state.ijk[K] == -1         && state.dir[Z] < 0.0)
           || (state.ijk[I] == extents[I] && state.dir[X] > 0.0)
           || (state.ijk[J] == extents[J] && state.dir[Y] > 0.0)
           || (state.ijk[K] == extents[K] && state.dir[Z] > 0.0))
    {
        state.next_dist = std::numeric_limits<double>::max();
        return state.next_dist;
    }

    const Vec_Dbl& edges_x(d_mesh.edges(X));
    const Vec_Dbl& edges_y(d_mesh.edges(Y));
    const Vec_Dbl& edges_z(d_mesh.edges(Z));

    // We want to transport to the furthest possible face, checking only faces
    // that are intersectable by the particle and closest to the particle
    double test_dist;

    // Check X faces if outside the extents
    if      (state.dir[X] < 0.0 && state.ijk[I] == extents[I])
    {
        test_dist = (edges_x.back() - state.pos[X]) / state.dir[X];
        if (test_dist > state.next_dist)
        {
            state.next_dist = test_dist;
            next_face       = Geo_State_t::PLUS_X;
        }
    }
    else if (state.dir[X] > 0.0 && state.ijk[I] == -1)
    {
        test_dist = (edges_x.front() - state.pos[X]) / state.dir[X];
        if (test_dist > state.next_dist)
        {
            state.next_dist = test_dist;
            next_face = Geo_State_t::MINUS_X;
        }
    }

    // Check Y faces if outside the extents
    if      (state.dir[Y] < 0.0 && state.ijk[J] == extents[J])
    {
        test_dist = (edges_y.back() - state.pos[Y]) / state.dir[Y];
        if (test_dist > state.next_dist)
        {
            state.next_dist = test_dist;
            next_face = Geo_State_t::PLUS_Y;
        }
    }
    else if (state.dir[Y] > 0.0 && state.ijk[J] == -1)
    {
        test_dist = (edges_y.front() - state.pos[Y]) / state.dir[Y];
        if (test_dist > state.next_dist)
        {
            state.next_dist = test_dist;
            next_face = Geo_State_t::MINUS_Y;
        }
    }

    // Check Z faces if outside the extents
    if      (state.dir[Z] < 0.0 && state.ijk[K] == extents[K])
    {
        test_dist = (edges_z.back() - state.pos[Z]) / state.dir[Z];
        if (test_dist > state.next_dist)
        {
            state.next_dist = test_dist;
            next_face = Geo_State_t::PLUS_Z;
        }
    }
    else if (state.dir[Z] > 0.0 && state.ijk[K] == -1)
    {
        test_dist = (edges_z.front() - state.pos[Z]) / state.dir[Z];
        if (test_dist > state.next_dist)
        {
            state.next_dist = test_dist;
            next_face = Geo_State_t::MINUS_Z;
        }
    }

    if (next_face == Geo_State_t::END_FACES)
    {
        // Particle is already in the interior
        state.next_ijk = state.ijk;
        state.next_dist = 0.;
    }
    else
    {
        // axes that we need to calculate new faces for
        size_t u, v;

        switch (next_face)
        {
            case Geo_State_t::MINUS_X:
                state.next_ijk[I] = 0;
                u = J;
                v = K;
                break;

            case Geo_State_t::PLUS_X:
                state.next_ijk[I] = extents[I] - 1;
                u = J;
                v = K;
                break;

            case Geo_State_t::MINUS_Y:
                state.next_ijk[J] = 0;
                u = I;
                v = K;
                break;

            case Geo_State_t::PLUS_Y:
                state.next_ijk[J] = extents[J] - 1;
                u = I;
                v = K;
                break;

            case Geo_State_t::MINUS_Z:
                state.next_ijk[K] = 0;
                u = J;
                v = I;
                break;

            case Geo_State_t::PLUS_Z:
                state.next_ijk[K] = extents[K] - 1;
                u = J;
                v = I;
                break;
            default:
                CHECK(0); // shouldn't ever get here
                break;
        }

        // Position after the next transport
        double nextu = state.pos[u] + state.next_dist * state.dir[u];
        double nextv = state.pos[v] + state.next_dist * state.dir[v];

        // Calculate IJK indices for the off-axis directions
        state.next_ijk[u] = d_mesh.find_upper(nextu, u);
        state.next_ijk[v] = d_mesh.find_upper(nextv, v);

        // Bump if on exact outer boundaries
        if (state.next_ijk[u] == -1 && soft_equiv(nextu, d_mesh.low_corner(u)))
            state.next_ijk[u] = 0;
        else if (state.next_ijk[u] == extents[u]
                 && soft_equiv(nextu, d_mesh.high_corner(u)))
            state.next_ijk[u] = extents[u] - 1;

        if (state.next_ijk[v] == -1 && soft_equiv(nextv, d_mesh.low_corner(v)))
            state.next_ijk[v] = 0;
        else if (state.next_ijk[v] == extents[v]
                 && soft_equiv(nextv, d_mesh.high_corner(v)))
            state.next_ijk[v] = extents[v] - 1;

        // If any of those are outside, we miss completely
        if (   state.next_ijk[u] < 0 || state.next_ijk[u] >= extents[u]
               || state.next_ijk[v] < 0 || state.next_ijk[v] >= extents[v])
        {
            state.next_dist = std::numeric_limits<double>::max();
        }
    }

    ENSURE(state.next_dist >= 0.);
    return state.next_dist;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create cell volumes and return them.
 *
 * On first pass, the volumes are generated; thereafter they are returned.
 */
Mesh_Geometry::SP_Vec_Dbl Mesh_Geometry::get_cell_volumes()
{
    using def::I; using def::J; using def::K;

    // if the volumes do not exist make them now (lazy evaluation)
    if (!d_volumes)
    {
        d_volumes     = std::make_shared<Vec_Dbl>(num_cells(), 0.0);
        auto &volumes = *d_volumes;
        const auto &N = d_mesh.extents();

        for (int k = 0; k < N[K]; ++k)
        {
            for (int j = 0; j < N[J]; ++j)
            {
                for (int i = 0; i < N[I]; ++i)
                {
                    volumes[d_mesh.index(i,j,k)] = d_mesh.volume(i,j,k);
                }
            }
        }
    }

    ENSURE(d_volumes->size() == num_cells());
    return d_volumes;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Reset tracking given current position, direction
 */
void Mesh_Geometry::update_state(Geo_State_t& state) const
{
    using def::I; using def::J; using def::K;

    // Find the logical indices of the cell along each axis
    d_mesh.find_upper(state.pos, state.ijk);

#ifdef CHECK_ON
    // Clear the other internals for debugging
    state.next_ijk[I] = -2;
    state.next_ijk[J] = -2;
    state.next_ijk[K] = -2;
    state.next_dist   = -1.;
#endif
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Mesh_Geometry.cc
//---------------------------------------------------------------------------//
