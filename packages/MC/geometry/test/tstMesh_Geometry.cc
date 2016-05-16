//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/geometry/test/tstMesh_Geometry.cc
 * \author Seth R Johnson
 * \date   Thu Mar 07 17:27:27 2013
 * \brief  Test for mesh geometry
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Mesh_Geometry.hh"

#include "utils/Definitions.hh"
#include "../Definitions.hh"

using profugus::Mesh_Geometry;
using def::I; using def::J; using def::K;
using def::X; using def::Y; using def::Z;
using profugus::geometry::INSIDE;
using profugus::geometry::OUTSIDE;
using profugus::geometry::REFLECT;

const double sqrt_third = 1./std::sqrt(3.0);
const double sqrt_half  = 1./std::sqrt(2.0);

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class MeshGeometryTest : public ::testing::Test
{
  protected:
    typedef Mesh_Geometry           Geometry_t;
    typedef Geometry_t::Geo_State_t Geo_State_t;
    typedef def::Vec_Dbl            Vec_Dbl;
    typedef def::Vec_Int            Vec_Int;

    typedef Geo_State_t::Space_Vector Space_Vector;
    typedef Geo_State_t::Coordinates  Coordinates;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        x = {0.0, 0.10, 0.25, 0.30, 0.42};
        y = {0.0, 0.20, 0.40, 0.50};
        z = {-0.1, 0.0, 0.15, 0.50};
    }

  protected:
    // >>> Data that get re-initialized between tests
    Vec_Dbl x;
    Vec_Dbl y;
    Vec_Dbl z;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, accessors)
{
    Geometry_t geo(x, y, z);

    EXPECT_EQ(36, geo.num_cells());

    // Test volume
    double expected_vol[] = {0.002, 0.003, 0.001, 0.0024, 0.002, 0.003, 0.001,
        0.0024, 0.001, 0.0015, 0.0005, 0.0012, 0.003, 0.0045, 0.0015, 0.0036,
        0.003, 0.0045, 0.0015, 0.0036, 0.0015, 0.00225, 0.00075, 0.0018, 0.007,
        0.0105, 0.0035, 0.0084, 0.007, 0.0105, 0.0035, 0.0084, 0.0035, 0.00525,
        0.00175, 0.0042};

    // Test all volumes
    auto volumes = geo.get_cell_volumes();
    ASSERT_TRUE(static_cast<bool>(volumes));

    EXPECT_VEC_SOFTEQ(expected_vol, *volumes, 1.e-12);
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, initialization)
{
    Mesh_Geometry geo(x, y, z);

    Geo_State_t state;

    // begin outside
    Space_Vector r(-1., 1., 0.1);
    Space_Vector omega(1., 0.0, 0.);
    geo.initialize(r, omega, state);

    // check cell index
#ifdef REQUIRE_ON
    EXPECT_THROW(geo.cell(state), profugus::assertion);
#endif
    EXPECT_EQ(OUTSIDE, geo.boundary_state(state));

    // check direction, position
    EXPECT_VEC_SOFT_EQ(r,     geo.position(state));
    EXPECT_VEC_SOFT_EQ(omega, geo.direction(state));

    // check state internals
    EXPECT_EQ(-1 , state.ijk[I]);
    EXPECT_EQ( 3 , state.ijk[J]);
    EXPECT_EQ( 1 , state.ijk[K]);

    // begin inside
    r = Space_Vector(0.2, 0.45, 0.4);
    omega = Space_Vector(0., 0., 1.);
    geo.initialize(r, omega, state);

    // check cell index
    EXPECT_EQ(1 + 2 * 4 + 2 * 4 * 3, geo.cell(state));
    Space_Vector pos = geo.position(state);
    EXPECT_EQ(1 + 2 * 4 + 2 * 4 * 3, geo.cell(pos));

    // check direction, position
    EXPECT_VEC_SOFT_EQ(r,     geo.position(state));
    EXPECT_VEC_SOFT_EQ(omega, geo.direction(state));

    // check state internals
    EXPECT_EQ(1, state.ijk[I]);
    EXPECT_EQ(2, state.ijk[J]);
    EXPECT_EQ(2, state.ijk[K]);
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, tags)
{
    Geometry_t geo(x, y, z);

    Geometry_t::SP_Vec_Int matids(std::make_shared<Vec_Int>(36));
    for (int i = 0; i < 36; ++i)
        (*matids)[i] = i * 2 + 1;

    // assign matids
    geo.set_matids(matids);

    Geo_State_t state;

    Space_Vector omega(1., 0.0, 0.);
    Space_Vector r;

    // Test tags
    int i,j,k;
    for (int cell = 0; cell < 36; ++cell)
    {
        geo.mesh().cardinal(cell, i, j, k);
        r[X] = 0.5 * (x[i] + x[i+1]);
        r[Y] = 0.5 * (y[j] + y[j+1]);
        r[Z] = 0.5 * (z[k] + z[k+1]);

        geo.initialize(r, omega, state);
        ASSERT_EQ(cell, geo.cell(state));
        EXPECT_EQ(cell        , geo.cell(state));
        EXPECT_EQ(cell * 2 + 1, geo.matid(state));

        EXPECT_EQ(INSIDE, geo.boundary_state(state));

        // Check matid search
        EXPECT_EQ(geo.matid(state), geo.matid(r));
    }

    r[Z] = 100.;
    geo.initialize(r, omega, state);
    EXPECT_EQ(OUTSIDE, geo.boundary_state(state));
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, distance_to_boundary)
{
    Mesh_Geometry geo(x, y, z);

    Geo_State_t state;
    Space_Vector r, omega;
    double dist;

    // distance from cell zero
    r = Space_Vector(0.01, 0.01, -0.01);
    omega = Space_Vector(1., 0., 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_boundary(state);
    EXPECT_SOFTEQ(0.1 - 0.01, dist, 1.e-12);
    EXPECT_EQ(Coordinates(1,0,0), state.next_ijk);

    // distance from another interior cell (slightly closer to X plane)
    r = Space_Vector(0.26, 0.35, -0.01);
    omega = Space_Vector(sqrt_half, sqrt_half, 0.);
    geo.initialize(r, omega, state);
    EXPECT_EQ(Coordinates(2,1,0), state.ijk);
    dist = geo.distance_to_boundary(state);
    EXPECT_SOFTEQ(0.04 / sqrt_half, dist, 1.e-12);
    EXPECT_EQ(Coordinates(3,1,0), state.next_ijk);

    // distance from cell zero, on boundary, initializes particle in cell -1
    // so we don't test that
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, move_to_point)
{
    Mesh_Geometry geo(x, y, z);

    Geo_State_t state;
    Space_Vector r, omega;

    // cell zero
    r = Space_Vector(0.01, 0.01, -0.01);
    omega = Space_Vector(1., 0., 0.);
    geo.initialize(r, omega, state);
    EXPECT_EQ(0, geo.cell(state));

    // move so it crosses cell boundaries
    geo.move_to_point(0.3, state);
    EXPECT_EQ(3, geo.cell(state));

    // move so it stays inside cell
    geo.move_to_point(0.00001, state);
    EXPECT_EQ(3, geo.cell(state));
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, distance_to_interior)
{
    Mesh_Geometry geo(x, y, z);

    Geo_State_t state;
    Space_Vector r, omega;
    double dist;

    // distance from the interior should be a null-op
    r = Space_Vector(0.01, 0.01, -0.01);
    omega = Space_Vector(1., 0., 0.);
    geo.initialize(r, omega, state);
    EXPECT_EQ(Coordinates(0,0,0), state.ijk);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(0., dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,0,0), state.next_ijk);

    // distance from exact -x boundary
    r = Space_Vector(0.0, 0.01, -0.01);
    omega = Space_Vector(1., 0., 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(0.0, dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,0,0), state.next_ijk);

    // lower corner case
    r = Space_Vector(-1.0, -1.0, -1.0 + -0.1);
    omega = Space_Vector(sqrt_third, sqrt_third, sqrt_third);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(std::sqrt(3.0), dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,0,0), state.next_ijk);

    // edge case
    r = Space_Vector(-1.0, -1.0, -1.0 + 0.01);
    omega = Space_Vector(sqrt_third, sqrt_third, sqrt_third);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(std::sqrt(3.0), dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,0,1), state.next_ijk);

    // upper corner case
    r = Space_Vector(1.42, 1.5, 1.5);
    omega = Space_Vector(-sqrt_third, -sqrt_third, -sqrt_third);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(std::sqrt(3.0), dist, 1.e-12);
    EXPECT_EQ(Coordinates(3,2,2), state.next_ijk);

    // distance from cell zero, left of boundary
    r = Space_Vector(-10.0, 0.01, -0.01);
    omega = Space_Vector(1., 0., 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(10., dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,0,0), state.next_ijk);

    // distance from positive X boundary
    r = Space_Vector(10.0, 0.01, -0.01);
    omega = Space_Vector(-1., 0., 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(10. - 0.42, dist, 1.e-12);
    EXPECT_EQ(Coordinates(3,0,0), state.next_ijk);

    // distance from negative Y boundary
    r = Space_Vector(0.35, -5., -0.01);
    omega = Space_Vector(0., 1., 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(5., dist, 1.e-12);
    EXPECT_EQ(Coordinates(3,0,0), state.next_ijk);

    // distance from positive Y boundary
    r = Space_Vector(0.01, 10.0, -0.01);
    omega = Space_Vector(0., -1., 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(10. - 0.5, dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,2,0), state.next_ijk);

    // distance from negative Z boundary
    r = Space_Vector(0.35, 0.01, -1.);
    omega = Space_Vector(0., 0., 1.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(1. - 0.1, dist, 1.e-12);
    EXPECT_EQ(Coordinates(3,0,0), state.next_ijk);

    // distance from positive Z boundary
    r = Space_Vector(0.01, 0.01, 100.);
    omega = Space_Vector(0., 0., -1.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_SOFTEQ(100. - 0.5, dist, 1.e-12);
    EXPECT_EQ(Coordinates(0,0,2), state.next_ijk);

    // misses the boundary by going over the cube
    r = Space_Vector(-2., -2., 2.);
    omega = Space_Vector(sqrt_third, sqrt_third, sqrt_third);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_GT(dist, 1.e10) << "new coordinates: " << state.next_ijk;

    // misses the boundary by z direction
    r = Space_Vector(0.01, 0.01, 100.);
    omega = Space_Vector(0., 0., 1.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_GT(dist, 1.e10) << "new coordinates: " << state.next_ijk;

    // misses the boundary in wrong direction
    r = Space_Vector(-2., -2., 0.);
    omega = Space_Vector(-sqrt_half, -sqrt_half, 0.);
    geo.initialize(r, omega, state);
    dist = geo.distance_to_interior(state);
    EXPECT_GT(dist, 1.e10);
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, tracking)
{
    using std::setw;

    Mesh_Geometry geo(x, y, z);

    Geo_State_t state;

    Space_Vector r(-0.01, 0.01, -0.01);
    //Space_Vector r(-10., -10., -10.);
    //Space_Vector r(10., -10., -10.);
    //Space_Vector omega(sqrt_half, sqrt_half, 0.);
    Space_Vector omega(sqrt_third, sqrt_third, sqrt_third);

    r      = Space_Vector( 0.4200,  0.2000,  0.0100);
    omega  = Space_Vector(-0.80178373, 0.53452248, 0.26726124);

    geo.initialize(r, omega, state);

    // Storage for cell index and distance traveled
    unsigned int cell = geo.cell(state);
    double dist;

    cout << "Initial condition:  "
            << "cell " << setw(4) << cell << " = " << state.ijk << "; "
            << "pos = "<< setw(10) << geo.position(state) << endl;

    // Calculate distance to edge of first cell
    dist = geo.distance_to_interior(state);
    //dist = 0.;

    if (dist > 0.)
    {
        // Move to inside edge and update cell
        geo.move_to_surface(state);
        cell = geo.cell(state);
        // If cell is -1, we've missed the tally mesh completely
        cout << "Moved " << setw(10) << dist << " to "
             << "cell " << setw(4) << cell << " = " << state.ijk << "; "
             << "pos = "<< setw(10) << geo.position(state) << endl;
    }

    while (geo.boundary_state(state) == INSIDE)
    {
        cell = geo.cell(state);
        // Calculate next distance
        dist = geo.distance_to_boundary(state);

        cout << "Moved " << setw(10) << dist << " in "
             << "cell " << setw(4) << cell << " = " << state.ijk << "; "
             << "pos = " << setw(10) << geo.position(state) << endl;

        // Move to next cell boundary and update cell
        geo.move_to_surface(state);
    }

    cout << "Done:" << "cell " << setw(4) << cell << " = " << state.ijk << "; "
         << "pos = " << geo.position(state) << endl;
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, cell_box)
{
    using def::I;
    using def::J;
    using def::K;

    Mesh_Geometry geo(x, y, z);

    int Nx = x.size() - 1;
    int Ny = y.size() - 1;
    int Nz = z.size() - 1;

    for( int k = 0; k < Nz; ++k )
    {
        for( int j = 0; j < Ny; ++j )
        {
            for( int i = 0; i < Nx; ++i )
            {
                int cell = i + Nx * (j + Ny * k);

                auto box = geo.get_cell_extents(cell);
                auto lower = box.lower();
                auto upper = box.upper();

                EXPECT_SOFT_EQ( lower[I], x[i] );
                EXPECT_SOFT_EQ( lower[J], y[j] );
                EXPECT_SOFT_EQ( lower[K], z[k] );
                EXPECT_SOFT_EQ( upper[I], x[i+1] );
                EXPECT_SOFT_EQ( upper[J], y[j+1] );
                EXPECT_SOFT_EQ( upper[K], z[k+1] );
            }
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(MeshGeometryTest, reflect)
{
    using def::I;
    using def::J;
    using def::K;

    Mesh_Geometry geo(x, y, z);

    std::vector<int> refl = {1, 0, 0, 0, 1, 1};
    geo.set_reflecting(refl);

    int Nx = x.size() - 1;
    int Ny = y.size() - 1;
    int Nz = z.size() - 1;

    //x = {0.0, 0.10, 0.25, 0.30, 0.42};
    //y = {0.0, 0.20, 0.40, 0.50};
    //z = {-0.1, 0.0, 0.15, 0.50};

    //
    // Test 1: reflection on reflecting surface
    //

    // Create state
    Geo_State_t state;
    def::Space_Vector r   = {0.05, 0.42, 0.10};
    def::Space_Vector dir = {-4.0, 0.1, -0.5};
    vector_normalize(dir);
    geo.initialize(r,dir,state);

    // Expect initial cell indices
    EXPECT_EQ( 0, state.ijk[I] );
    EXPECT_EQ( 2, state.ijk[J] );
    EXPECT_EQ( 1, state.ijk[K] );

    // Move to surface
    geo.distance_to_boundary(state);
    geo.move_to_surface(state);

    EXPECT_EQ( Geo_State_t::MINUS_X, state.exiting_face );
    EXPECT_EQ( Geo_State_t::MINUS_X, state.reflecting_face);
    EXPECT_EQ( REFLECT,              geo.boundary_state(state) );

    // Cell indices should be unchanged
    EXPECT_EQ( 0, state.ijk[I] );
    EXPECT_EQ( 2, state.ijk[J] );
    EXPECT_EQ( 1, state.ijk[K] );

    // Now reflect
    bool reflected = geo.reflect(state);
    EXPECT_TRUE( reflected );

    // Check direction
    dir[I] = -dir[I];
    EXPECT_SOFT_EQ( dir[I], state.d_dir[I] );
    EXPECT_SOFT_EQ( dir[J], state.d_dir[J] );
    EXPECT_SOFT_EQ( dir[K], state.d_dir[K] );

    // Cell index shouldn't change during reflection
    EXPECT_EQ( 0, state.ijk[I] );
    EXPECT_EQ( 2, state.ijk[J] );
    EXPECT_EQ( 1, state.ijk[K] );

    //
    // Test 2: reflection on non-reflecting surface
    //

    // Create state
    r   = {0.20, 0.45, 0.35};
    dir = {-1.0, 2.0, 1.0};
    vector_normalize(dir);
    geo.initialize(r,dir,state);

    // Initial cell indices
    EXPECT_EQ( 1,  state.ijk[I] );
    EXPECT_EQ( 2, state.ijk[J] );
    EXPECT_EQ( 2,  state.ijk[K] );

    // Move to surface
    geo.distance_to_boundary(state);
    geo.move_to_surface(state);

    EXPECT_EQ( Geo_State_t::PLUS_Y,  state.exiting_face );
    EXPECT_EQ( Geo_State_t::NONE,    state.reflecting_face);
    EXPECT_EQ( OUTSIDE,              geo.boundary_state(state) );

    // Y index should be changed
    EXPECT_EQ( 1,  state.ijk[I] );
    EXPECT_EQ( Ny, state.ijk[J] );
    EXPECT_EQ( 2,  state.ijk[K] );

    // Test reflection on non-reflecting face
    reflected = geo.reflect(state);
    EXPECT_FALSE( reflected );

    // Check direction
    EXPECT_SOFT_EQ( dir[I], state.d_dir[I] );
    EXPECT_SOFT_EQ( dir[J], state.d_dir[J] );
    EXPECT_SOFT_EQ( dir[K], state.d_dir[K] );

    EXPECT_EQ( 1,  state.ijk[I] );
    EXPECT_EQ( Ny, state.ijk[J] );
    EXPECT_EQ( 2,  state.ijk[K] );

}

//---------------------------------------------------------------------------//
//                        end of tstMesh_Geometry.cc
//---------------------------------------------------------------------------//
