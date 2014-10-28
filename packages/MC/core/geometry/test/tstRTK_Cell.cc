//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/test/tstRTK_Cell.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 21 13:16:16 2010
 * \brief  RTK_Cell unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "utils/Constants.hh"
#include "rng/RNG_Control.hh"
#include "../RTK_Cell.hh"

using namespace std;

using profugus::RTK_Cell;
using profugus::RNG_Control;

typedef RTK_Cell::Space_Vector Vector;
typedef RTK_Cell::Geo_State_t  Geo_State;

int seed = 1235123;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(State, Packing)
{
    // make a buffer
    vector<char> buffer;

    EXPECT_TRUE(Geo_State::packed_bytes() == 16 * sizeof(int) +
              6 * sizeof(double));

    // pack a state
    {
        Geo_State state;

        state.d_r[0]   = 1.1;
        state.d_r[1]   = 1.2;
        state.d_r[2]   = 1.3;
        state.d_dir[0] = 0.4;
        state.d_dir[1] = 0.1;
        state.d_dir[2] = 0.2;

        state.region       = 4;
        state.segment      = 2;
        state.face         = Geo_State::NONE;
        state.next_face    = 1;
        state.next_region  = 2;
        state.next_segment = 3;

        state.escaping_face = Geo_State::NONE;
        state.exiting_face  = Geo_State::MINUS_Y;

        state.level_coord[0][0] = 1;
        state.level_coord[0][1] = 2;
        state.level_coord[0][2] = 3;

        state.level_coord[1][0] = 4;
        state.level_coord[1][1] = 5;
        state.level_coord[1][2] = 6;

        state.level_coord[2][0] = 7;
        state.level_coord[2][1] = 8;
        state.level_coord[2][2] = 9;

        // assume 4 bits of junk at the beginning
        buffer.resize(4 + state.packed_bytes());

        // pack the data
        state.pack(&buffer[4]);
    }

    // unpack the state
    {
        Geo_State state;
        state.unpack(&buffer[4]);

        EXPECT_EQ(1.1, state.d_r[0]);
        EXPECT_EQ(1.2, state.d_r[1]);
        EXPECT_EQ(1.3, state.d_r[2]);
        EXPECT_EQ(0.4, state.d_dir[0]);
        EXPECT_EQ(0.1, state.d_dir[1]);
        EXPECT_EQ(0.2, state.d_dir[2]);

        EXPECT_EQ(4, state.region);
        EXPECT_EQ(Geo_State::NONE, state.face);
        EXPECT_EQ(2, state.segment);
        EXPECT_EQ(1, state.next_face);
        EXPECT_EQ(2, state.next_region);
        EXPECT_EQ(3, state.next_segment);
        EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);

        EXPECT_EQ(1, state.level_coord[0][0]);
        EXPECT_EQ(2, state.level_coord[0][1]);
        EXPECT_EQ(3, state.level_coord[0][2]);

        EXPECT_EQ(4, state.level_coord[1][0]);
        EXPECT_EQ(5, state.level_coord[1][1]);
        EXPECT_EQ(6, state.level_coord[1][2]);

        EXPECT_EQ(7, state.level_coord[2][0]);
        EXPECT_EQ(8, state.level_coord[2][1]);
        EXPECT_EQ(9, state.level_coord[2][2]);
    }
}

//---------------------------------------------------------------------------//

TEST(Single, Shell)
{
    // make a pin cell
    RTK_Cell pin1(1, 0.54, 10, 1.26, 14.28);
    RTK_Cell pin2(2, 0.46, 10, 1.26, 14.28);
    RTK_Cell pin4(4, 0.63, 10, 1.26, 14.28);

    EXPECT_EQ(2, pin1.num_regions());
    EXPECT_EQ(1, pin1.num_shells());
    EXPECT_EQ(2, pin2.num_regions());
    EXPECT_EQ(1, pin2.num_shells());
    EXPECT_EQ(2, pin4.num_regions());
    EXPECT_EQ(1, pin4.num_shells());

    using def::X; using def::Y; using def::Z;
    Vector lower, upper;
    pin1.get_extents(lower, upper);

    EXPECT_SOFTEQ(lower[X], -1.26/2, 1.e-12);
    EXPECT_SOFTEQ(lower[Y], -1.26/2, 1.e-12);
    EXPECT_SOFTEQ(lower[Z],  0.    , 1.e-12);
    EXPECT_SOFTEQ(upper[X],  1.26/2, 1.e-12);
    EXPECT_SOFTEQ(upper[Y],  1.26/2, 1.e-12);
    EXPECT_SOFTEQ(upper[Z], 14.28  , 1.e-12);

    // test some queries to find the shell a point is in
    EXPECT_EQ(1, pin1.region(0.0, 0.5401));
    EXPECT_EQ(0, pin1.region(0.0, 0.5399));
    EXPECT_EQ(1, pin2.region(0.4601, 0.0));
    EXPECT_EQ(0, pin2.region(0.4600, 0.0));
    EXPECT_EQ(0, pin2.region(0.4599, 0.0));
    EXPECT_EQ(1, pin4.region(-0.62, 0.12));
    EXPECT_EQ(0, pin4.region(-0.62, 0.11));

    EXPECT_EQ(1, pin1.matid(0));
    EXPECT_EQ(10, pin1.matid(1));
    EXPECT_EQ(2, pin2.matid(0));
    EXPECT_EQ(10, pin2.matid(1));
    EXPECT_EQ(4, pin4.matid(0));
    EXPECT_EQ(10, pin4.matid(1));

#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        EXPECT_TRUE(pin4.region(.6301, 0.0));
    }
    catch(const profugus::assertion &ass)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif

    // test initialization
    Geo_State state;
    pin1.initialize(Vector(0.0, 0.53, 0.0), state);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(0, state.segment);
    pin1.initialize(Vector(0.0, 0.54, 0.0), state);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(0, state.segment);
    EXPECT_EQ(0, pin1.cell(state.region, state.segment));
    pin1.initialize(Vector(0.0, 0.55, 0.0), state);
    EXPECT_EQ(1, state.region);
    EXPECT_EQ(0, state.segment);
    EXPECT_EQ(1, pin1.cell(state.region, state.segment));

    // distance to boundary

    // box-boundary tests
    Vector r, omega;
    {
        r     = Vector(0.0, 0.59, 0.0);
        omega = Vector(1.0, 0.0, 0.0);
        pin1.initialize(r, state);
        pin1.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.63, 1.e-12);
        EXPECT_EQ(Geo_State::PLUS_X, state.exiting_face);
        EXPECT_EQ(1, state.region);
    }
    {
        r     = Vector(0.0, 0.59, 0.0);
        omega = Vector(-1.0, 0.0, 0.0);
        pin1.initialize(r, state);
        pin1.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.63, 1.e-12);
        EXPECT_EQ(Geo_State::MINUS_X, state.exiting_face);
        EXPECT_EQ(1, state.region);
    }

    // Pin intersections
    RTK_Cell pin(1, 0.45, 2, 1.2, 14.28);
    double   eps = 1.0e-6;
    {
        r     = Vector(  0.43,   0.51,   1.20);
        omega = Vector( -0.07450781,  -0.17272265,   0.98214840);
        pin.initialize(r, state);
        EXPECT_EQ(1, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.2334036420, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);
    }
    {
        r     = Vector(  0.43,   0.51,   1.20);
        omega = Vector(  0.01923789,  -0.98113214,  -0.19237885);
        pin.initialize(r, state);
        EXPECT_EQ(1, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.41448110826, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);
    }
    {
        r     = Vector( -0.49,   0.10,   1.20);
        omega = Vector(  0.04377546,  -0.01122448,   0.99897834);
        pin.initialize(r, state);
        EXPECT_EQ(1, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.1101358881, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);
    }
    {
        r     = Vector( -0.31,  -0.46,   1.20);
        omega = Vector( -0.01642565,   0.05812155,   0.99817438);
        pin.initialize(r, state);
        EXPECT_EQ(1, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 3.4103552300, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);
    }
    {
        r     = Vector(  0.21,  -0.56,   1.20);
        omega = Vector( -0.03911234,   0.07065454,   0.99673375);
        pin.initialize(r, state);
        EXPECT_EQ(1, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.860292469, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);
    }
    {
        r     = Vector(  0.21,  -0.56,  14.20);
        omega = Vector(  0.021262916894,   0.051770580263,   0.998432619351);
        pin.initialize(r, state);
        EXPECT_EQ(1, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.08 / omega[2], eps);
        EXPECT_EQ(Geo_State::PLUS_Z, state.exiting_face);
    }

    {
        r     = Vector(  0.42,   0.10,  14.20);
        omega = Vector( -0.262232636986,  -0.030141682412,  -0.964533837188);
        pin.initialize(r, state);
        EXPECT_EQ(0, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 3.3176648414, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.next_region);
    }
    {
        r     = Vector(  0.10,   0.30,  12.10);
        omega = Vector(  0.408248290464,   0.163299316186,  -0.898146239020);
        pin.initialize(r, state);
        EXPECT_EQ(0, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 3.9914694397e-01, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.next_region);
    }
    {
        r     = Vector(  0.10,   0.30,  12.10);
        omega = Vector(  0.011875513070,   0.004750205228,  -0.999918200524);
        pin.initialize(r, state);
        EXPECT_EQ(0, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, -12.10 / omega[2], eps);
        EXPECT_EQ(Geo_State::MINUS_Z, state.exiting_face);
    }
    {
        r     = Vector(  0.10,   0.30,  12.10);
        omega = Vector(  0.072244077132,   0.028897630853,   0.996968264415);
        pin.initialize(r, state);
        EXPECT_EQ(0, state.region);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 2.18 / omega[2], eps);
        EXPECT_EQ(Geo_State::PLUS_Z, state.exiting_face);
    }
}

//---------------------------------------------------------------------------//

TEST(Multi, Shell)
{
    // make pin with clad
    vector<int>    ids(2, 0);
    vector<double> rad(2, 0);
    ids[0] = 1;
    rad[0] = 0.49;
    ids[1] = 2;
    rad[1] = 0.54;

    // Pin
    RTK_Cell pin(ids, rad, 3, 1.26, 14.28);

    EXPECT_SOFT_EQ(0.49, pin.radii()[0]);
    EXPECT_SOFT_EQ(0.54, pin.radii()[1]);

    EXPECT_EQ(3, pin.num_regions());
    EXPECT_EQ(2, pin.num_shells());
    EXPECT_EQ(1, pin.num_segments());
    EXPECT_EQ(1.26, pin.pitch(0));
    EXPECT_EQ(1.26, pin.pitch(1));
    EXPECT_EQ(14.28, pin.height());

    EXPECT_EQ(1, pin.matid(0));
    EXPECT_EQ(2, pin.matid(1));
    EXPECT_EQ(3, pin.matid(2));

    EXPECT_EQ(1, pin.region(0.1, 0.48));
    EXPECT_EQ(0, pin.region(0.1, 0.479));
    EXPECT_EQ(2, pin.region(0.5, 0.3));
    EXPECT_EQ(1, pin.region(0.4, 0.35));

    // tracking test
    Geo_State state;
    pin.initialize(Vector(0.5, 0.3, 12.1), state);
    EXPECT_EQ(2, state.region);
    EXPECT_EQ(Geo_State::NONE, state.face);

    Vector r, omega;
    double eps = 1.0e-6;
    {
        r     = Vector(  0.50,   0.30,  12.10);
        omega = Vector( -0.740797197487,  -0.642024237822,   0.197545919330);
        pin.initialize(r, state);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.0446878772402, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(2, state.region);
        EXPECT_EQ(1, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.0520128055639, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(0, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.978336739656, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.region);
        EXPECT_EQ(1, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.0520128055639, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(2, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(2, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.32149321506, eps);
        EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);
        EXPECT_EQ(2, state.region);
    }
}

//---------------------------------------------------------------------------//

TEST(Empty, SquareCell)
{
    // make an empty pin cell
    RTK_Cell box(11, 1.26, 14.28);

    EXPECT_EQ(1, box.num_regions());
    EXPECT_EQ(0, box.num_shells());

    EXPECT_EQ(0, box.region(0.0, 0.5401));
    EXPECT_EQ(0, box.region(0.0, 0.5399));

    EXPECT_EQ(11, box.matid(0));

    Geo_State state;
    box.initialize(Vector(0.0, 0.53, 0.0), state);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(Geo_State::NONE, state.face);

    Vector r, omega;
    double eps = 1.0e-6;
    {
        r     = Vector(0.0, 0.59, 0.0);
        omega = Vector(1.0, 0.0, 0.0);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.63, 1.e-12);
        EXPECT_EQ(Geo_State::PLUS_X, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector(0.0, 0.59, 0.0);
        omega = Vector(-1.0, 0.0, 0.0);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.63, 1.e-12);
        EXPECT_EQ(Geo_State::MINUS_X, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector(  0.23,  -0.63,  12.10);
        omega = Vector( -0.591113929288,   0.783346101414,   0.192232172126);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.4548802818, eps);
        EXPECT_EQ(Geo_State::MINUS_X, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector(  0.45,  -0.62,  12.10);
        omega = Vector( -0.628969195431,   0.754763034517,   0.186361243091);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.6561489406, eps);
        EXPECT_EQ(Geo_State::PLUS_Y, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector( -0.12,  -0.38,  12.10);
        omega = Vector(  0.810238620974,  -0.492497985298,   0.317740635676);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.50761628974, eps);
        EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector( -0.12,  -0.38,  12.10);
        omega = Vector(  0.097437248619,  -0.059226562886,   0.993477829058);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 2.1943116758, eps);
        EXPECT_EQ(Geo_State::PLUS_Z, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector( -0.12,  -0.38,   2.10);
        omega = Vector(  0.080591552144,  -0.048987021891,  -0.995542702956);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 2.1094022323, eps);
        EXPECT_EQ(Geo_State::MINUS_Z, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
}

//---------------------------------------------------------------------------//

TEST(Empty, Cell)
{
    // make an empty pin cell
    RTK_Cell box(11, 0.25, 1.26, 14.28);

    EXPECT_EQ(1, box.num_regions());
    EXPECT_EQ(0, box.num_shells());

    EXPECT_EQ(0, box.region(0.0, 0.5401));
    EXPECT_EQ(0, box.region(0.0, 0.5399));

    EXPECT_EQ(11, box.matid(0));

    EXPECT_EQ(0.25, box.pitch(0));
    EXPECT_EQ(1.26, box.pitch(1));

    Geo_State state;
    box.initialize(Vector(0.0, 0.53, 0.0), state);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(Geo_State::NONE, state.face);

    Vector r, omega;
    double eps = 1.0e-6;
    {
        r     = Vector(0.0, 0.59, 0.0);
        omega = Vector(1.0, 0.0, 0.0);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.125, 1.e-12);
        EXPECT_EQ(Geo_State::PLUS_X, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector(0.0, 0.59, 0.0);
        omega = Vector(-1.0, 0.0, 0.0);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.125, 1.e-12);
        EXPECT_EQ(Geo_State::MINUS_X, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector(  0.123,  -0.63,  12.10);
        omega = Vector( -0.556102184819,   0.807165237093,   0.198077358796);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 4.4596120420e-01, eps);
        EXPECT_EQ(Geo_State::MINUS_X, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
    {
        r     = Vector(  0.045,  -0.62,  12.10);
        omega = Vector( -0.080640139920,   0.967681679037,   0.238933747910);
        box.initialize(r, state);
        box.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.2917470973e+00, eps);
        EXPECT_EQ(Geo_State::PLUS_Y, state.exiting_face);
        EXPECT_EQ(0, state.region);
    }
}

//---------------------------------------------------------------------------//
// See support/pin_cell_test.png for pin-cell face/cell mapping

TEST(Multi, SegPin)
{
    // make pin with clad and 4 segments
    vector<int>    ids(3, 0);
    vector<double> rad(3, 0);
    ids[0] = 1;
    rad[0] = 0.27;
    ids[1] = 2;
    rad[1] = 0.486;
    ids[2] = 5;
    rad[2] = 0.54;

    // Pin
    RTK_Cell pin(ids, rad, 10, 1.26, 14.28, 4);

    EXPECT_EQ(16, pin.num_cells());
    EXPECT_EQ(4, pin.num_regions());
    EXPECT_EQ(3, pin.num_shells());
    EXPECT_EQ(4, pin.num_segments());

    EXPECT_EQ(1, pin.region(0.28, 0.0));
    EXPECT_EQ(0, pin.segment(0.28, 0.0));
    EXPECT_EQ(2, pin.segment(0.28, -0.01));

    EXPECT_EQ(0, pin.cell(0, 0));
    EXPECT_EQ(1, pin.cell(1, 0));
    EXPECT_EQ(2, pin.cell(2, 0));
    EXPECT_EQ(3, pin.cell(3, 0));

    EXPECT_EQ(4, pin.cell(0, 1));
    EXPECT_EQ(5, pin.cell(1, 1));
    EXPECT_EQ(6, pin.cell(2, 1));
    EXPECT_EQ(7, pin.cell(3, 1));

    EXPECT_EQ(8, pin.cell(0, 2));
    EXPECT_EQ(9, pin.cell(1, 2));
    EXPECT_EQ(10, pin.cell(2, 2));
    EXPECT_EQ(11, pin.cell(3, 2));

    EXPECT_EQ(12, pin.cell(0, 3));
    EXPECT_EQ(13, pin.cell(1, 3));
    EXPECT_EQ(14, pin.cell(2, 3));
    EXPECT_EQ(15, pin.cell(3, 3));

    Vector r, omega;
    double eps = 1.0e-6;

    Geo_State state;

    // test tracking
    {
        r     = Vector(  0.43,   0.51,   1.20);
        omega = Vector( -0.267261241912,  -0.534522483825,   0.801783725737);

        // initialize
        pin.initialize(r, state);
        EXPECT_EQ(3, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(Geo_State::NONE, state.face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 2.2028008712e-01, eps);
        EXPECT_EQ(Geo_State::NONE, state.face);
        EXPECT_EQ(3, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(2, state.next_face);
        EXPECT_EQ(2, state.next_region);
        EXPECT_EQ(0, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(2, state.face);
        EXPECT_EQ(2, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 9.4898743120e-02, eps);
        EXPECT_EQ(2, state.face);
        EXPECT_EQ(2, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(1, state.next_face);
        EXPECT_EQ(1, state.next_region);
        EXPECT_EQ(0, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(1, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 4.0177140025e-01, eps);
        EXPECT_EQ(1, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(0, state.next_face);
        EXPECT_EQ(0, state.next_region);
        EXPECT_EQ(0, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(0, state.face);
        EXPECT_EQ(0, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 2.3717240314e-01, eps);
        EXPECT_EQ(0, state.face);
        EXPECT_EQ(0, state.region);
        EXPECT_EQ(0, state.segment);
        EXPECT_EQ(4, state.next_face);
        EXPECT_EQ(0, state.next_region);
        EXPECT_EQ(2, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(4, state.face);
        EXPECT_EQ(0, state.region);
        EXPECT_EQ(2, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 4.9908842021e-01, eps);
        EXPECT_EQ(4, state.face);
        EXPECT_EQ(0, state.region);
        EXPECT_EQ(2, state.segment);
        EXPECT_EQ(0, state.next_face);
        EXPECT_EQ(1, state.next_region);
        EXPECT_EQ(2, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(0, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(2, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.5570162247e-01, eps);
        EXPECT_EQ(0, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(2, state.segment);
        EXPECT_EQ(3, state.next_face);
        EXPECT_EQ(1, state.next_region);
        EXPECT_EQ(3, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(3, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(3, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 2.4606977777e-01, eps);
        EXPECT_EQ(3, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(3, state.segment);
        EXPECT_EQ(1, state.next_face);
        EXPECT_EQ(2, state.next_region);
        EXPECT_EQ(3, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(1, state.face);
        EXPECT_EQ(2, state.region);
        EXPECT_EQ(3, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 9.4898743120e-02, eps);
        EXPECT_EQ(1, state.face);
        EXPECT_EQ(2, state.region);
        EXPECT_EQ(3, state.segment);
        EXPECT_EQ(2, state.next_face);
        EXPECT_EQ(3, state.next_region);
        EXPECT_EQ(3, state.next_segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // move to surface
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        // cross the surface
        pin.cross_surface(state);

        EXPECT_EQ(2, state.face);
        EXPECT_EQ(3, state.region);
        EXPECT_EQ(3, state.segment);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

        // get distance to boundary
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 1.8286351326e-01, eps);
        EXPECT_EQ(2, state.face);
        EXPECT_EQ(3, state.region);
        EXPECT_EQ(3, state.segment);
        EXPECT_EQ(Geo_State::NONE, state.next_face);
        EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);
    }

    RNG_Control control(seed);
    auto rng = control.rng();

    // sample some points in the pin
    int N           = 1000;
    double costheta = 0.0;
    double sintheta = 0.0;
    double phi      = 0.0;
    for (int n = 0; n < N; ++n)
    {
        r[0] = 1.26 * rng.ran() - 0.63;
        r[1] = 1.26 * rng.ran() - 0.63;
        r[2] = 14.28 * rng.ran();

        costheta = 1.0 - 2.0 * rng.ran();
        sintheta = sqrt(1.0 - costheta * costheta);
        phi      = 2.0 * profugus::constants::pi * rng.ran();

        omega[0] = sintheta * cos(phi);
        omega[1] = sintheta * sin(phi);
        omega[2] = costheta;

        pin.initialize(r, state);

        EXPECT_TRUE(state.region >= 0 && state.region < 4);

        // track distance to boundary
        while (state.exiting_face == Geo_State::INTERNAL)
        {
            pin.distance_to_boundary(r, omega, state);

            // move to surface
            r[0] += omega[0] * state.dist_to_next_region;
            r[1] += omega[1] * state.dist_to_next_region;
            r[2] += omega[2] * state.dist_to_next_region;

            // cross the surface
            pin.cross_surface(state);
        }
    }
}

//---------------------------------------------------------------------------//

TEST(Gap, LoX_LoY)
{
    // make pin with clad
    vector<int>    ids(2, 0);
    vector<double> rad(2, 0);
    ids[0] = 1;
    rad[0] = 0.49;
    ids[1] = 2;
    rad[1] = 0.54;

    // add gap to low xy
    RTK_Cell::Gap_Vector gap(0.1, 0.0, 0.1, 0.0);

    // Pin
    RTK_Cell pin(ids, rad, 3, 1.26, 14.28, gap);

    Vector l, u;
    pin.get_extents(l, u);
    EXPECT_SOFT_EQ(-0.73, l[0]);
    EXPECT_SOFT_EQ(-0.73, l[1]);
    EXPECT_SOFT_EQ(0.0,   l[2]);
    EXPECT_SOFT_EQ(0.63,  u[0]);
    EXPECT_SOFT_EQ(0.63,  u[1]);
    EXPECT_SOFT_EQ(14.28, u[2]);

    EXPECT_SOFT_EQ(0.49, pin.radii()[0]);
    EXPECT_SOFT_EQ(0.54, pin.radii()[1]);

    EXPECT_EQ(3, pin.num_regions());
    EXPECT_EQ(2, pin.num_shells());
    EXPECT_EQ(1, pin.num_segments());
    EXPECT_SOFT_EQ(1.36, pin.pitch(0));
    EXPECT_SOFT_EQ(1.36, pin.pitch(1));
    EXPECT_EQ(14.28, pin.height());

    EXPECT_EQ(1, pin.matid(0));
    EXPECT_EQ(2, pin.matid(1));
    EXPECT_EQ(3, pin.matid(2));

    EXPECT_EQ(1, pin.region(0.1, 0.48));
    EXPECT_EQ(0, pin.region(0.1, 0.479));
    EXPECT_EQ(2, pin.region(0.5, 0.3));
    EXPECT_EQ(1, pin.region(0.4, 0.35));
    EXPECT_EQ(2, pin.region(-0.7, -0.7));
    EXPECT_EQ(2, pin.region(-0.64, 0.1));
    EXPECT_EQ(2, pin.region(0.1, -0.64));

#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        pin.region(.6301, 0.0);
    }
    catch(const profugus::assertion &ass)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif

    // tracking test
    Geo_State state;
    pin.initialize(Vector(0.5, 0.3, 12.1), state);
    EXPECT_EQ(2, state.region);
    EXPECT_EQ(Geo_State::NONE, state.face);

    Vector r, omega;
    double eps = 1.0e-6;
    {
        r     = Vector(  0.50,   0.30,  12.10);
        omega = Vector( -0.740797197487,  -0.642024237822,   0.197545919330);
        pin.initialize(r, state);
        pin.distance_to_boundary(r, omega, state);
        EXPECT_SOFTEQ(state.dist_to_next_region, 0.0446878772402, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(2, state.region);
        EXPECT_EQ(1, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.0520128055639, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(0, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.978336739656, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.region);
        EXPECT_EQ(1, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.0520128055639, eps);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(2, state.next_region);

        pin.cross_surface(state);
        EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
        EXPECT_EQ(2, state.region);

        for (int i = 0; i < 3; ++i)
            r[i] += state.dist_to_next_region * omega[i];

        pin.distance_to_boundary(r, omega, state);

        EXPECT_SOFTEQ(state.dist_to_next_region, 0.4772505785762855, eps);
        EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);
        EXPECT_EQ(2, state.region);
    }
}

//---------------------------------------------------------------------------//

TEST(Gap, LoX_HiY)
{
    // make pin with clad
    vector<int>    ids(2, 0);
    vector<double> rad(2, 0);
    ids[0] = 1;
    rad[0] = 0.49;
    ids[1] = 2;
    rad[1] = 0.54;

    // add gap to low xy
    RTK_Cell::Gap_Vector gap(0.1, 0.0, 0.0, 0.1);

    // Pin
    RTK_Cell pin(ids, rad, 3, 1.26, 14.28, gap);

    Vector l, u;
    pin.get_extents(l, u);
    EXPECT_SOFT_EQ(-0.73, l[0]);
    EXPECT_SOFT_EQ(-0.63, l[1]);
    EXPECT_SOFT_EQ(0.0,   l[2]);
    EXPECT_SOFT_EQ(0.63,  u[0]);
    EXPECT_SOFT_EQ(0.73,  u[1]);
    EXPECT_SOFT_EQ(14.28, u[2]);

    EXPECT_SOFT_EQ(0.49, pin.radii()[0]);
    EXPECT_SOFT_EQ(0.54, pin.radii()[1]);

    EXPECT_EQ(3, pin.num_regions());
    EXPECT_EQ(2, pin.num_shells());
    EXPECT_EQ(1, pin.num_segments());
    EXPECT_SOFT_EQ(1.36, pin.pitch(0));
    EXPECT_SOFT_EQ(1.36, pin.pitch(1));
    EXPECT_EQ(14.28, pin.height());

    EXPECT_EQ(1, pin.matid(0));
    EXPECT_EQ(2, pin.matid(1));
    EXPECT_EQ(3, pin.matid(2));

    EXPECT_EQ(1, pin.region(0.1, 0.48));
    EXPECT_EQ(0, pin.region(0.1, 0.479));
    EXPECT_EQ(2, pin.region(0.5, 0.3));
    EXPECT_EQ(1, pin.region(0.4, 0.35));
    EXPECT_EQ(2, pin.region(-0.7, 0.7));
    EXPECT_EQ(2, pin.region(-0.64, 0.1));
    EXPECT_EQ(2, pin.region(0.1, 0.64));

#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        pin.region(0.0, -0.6301);
    }
    catch(const profugus::assertion &ass)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}

//---------------------------------------------------------------------------//

TEST(Gap, HiX_HiY)
{
    // make pin with clad
    vector<int>    ids(2, 0);
    vector<double> rad(2, 0);
    ids[0] = 1;
    rad[0] = 0.49;
    ids[1] = 2;
    rad[1] = 0.54;

    // add gap to low xy
    RTK_Cell::Gap_Vector gap(0.0, 0.1, 0.0, 0.1);

    // Pin
    RTK_Cell pin(ids, rad, 3, 1.26, 14.28, gap);

    Vector l, u;
    pin.get_extents(l, u);
    EXPECT_SOFT_EQ(-0.63, l[0]);
    EXPECT_SOFT_EQ(-0.63, l[1]);
    EXPECT_SOFT_EQ(0.0,   l[2]);
    EXPECT_SOFT_EQ(0.73,  u[0]);
    EXPECT_SOFT_EQ(0.73,  u[1]);
    EXPECT_SOFT_EQ(14.28, u[2]);

    EXPECT_SOFT_EQ(0.49, pin.radii()[0]);
    EXPECT_SOFT_EQ(0.54, pin.radii()[1]);

    EXPECT_EQ(3, pin.num_regions());
    EXPECT_EQ(2, pin.num_shells());
    EXPECT_EQ(1, pin.num_segments());
    EXPECT_SOFT_EQ(1.36, pin.pitch(0));
    EXPECT_SOFT_EQ(1.36, pin.pitch(1));
    EXPECT_EQ(14.28, pin.height());

    EXPECT_EQ(1, pin.matid(0));
    EXPECT_EQ(2, pin.matid(1));
    EXPECT_EQ(3, pin.matid(2));

    EXPECT_EQ(1, pin.region(0.1, 0.48));
    EXPECT_EQ(0, pin.region(0.1, 0.479));
    EXPECT_EQ(2, pin.region(0.5, 0.3));
    EXPECT_EQ(1, pin.region(0.4, 0.35));
    EXPECT_EQ(2, pin.region(0.7, 0.7));
    EXPECT_EQ(2, pin.region(0.64, 0.1));
    EXPECT_EQ(2, pin.region(0.1, 0.64));

#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        pin.region(-0.6301, -0.6301);
    }
    catch(const profugus::assertion &ass)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}

//---------------------------------------------------------------------------//

TEST(GapHomCell, HiX_HiY)
{
    // add gap to low xy
    RTK_Cell::Gap_Vector gap(0.0, 0.1, 0.0, 0.1);

    // Pin
    RTK_Cell pin(3, 1.26, 14.28, gap);

    Vector l, u;
    pin.get_extents(l, u);
    EXPECT_SOFT_EQ(-0.63, l[0]);
    EXPECT_SOFT_EQ(-0.63, l[1]);
    EXPECT_SOFT_EQ(0.0,   l[2]);
    EXPECT_SOFT_EQ(0.73,  u[0]);
    EXPECT_SOFT_EQ(0.73,  u[1]);
    EXPECT_SOFT_EQ(14.28, u[2]);

    EXPECT_EQ(1, pin.num_regions());
    EXPECT_EQ(0, pin.num_shells());
    EXPECT_EQ(1, pin.num_segments());
    EXPECT_SOFT_EQ(1.36, pin.pitch(0));
    EXPECT_SOFT_EQ(1.36, pin.pitch(1));
    EXPECT_EQ(14.28, pin.height());

    EXPECT_EQ(3, pin.matid(0));

    EXPECT_EQ(0, pin.region(0.1, 0.48));
    EXPECT_EQ(0, pin.region(0.1, 0.479));
    EXPECT_EQ(0, pin.region(0.5, 0.3));
    EXPECT_EQ(0, pin.region(0.4, 0.35));
    EXPECT_EQ(0, pin.region(0.7, 0.7));
    EXPECT_EQ(0, pin.region(0.64, 0.1));
    EXPECT_EQ(0, pin.region(0.1, 0.64));

#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        pin.region(-0.6301, -0.6301);
    }
    catch(const profugus::assertion &ass)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}

//---------------------------------------------------------------------------//

TEST(Vessel, Add)
{
    // make a pins
    RTK_Cell water(10, 21.62, 30.0, 14.28);
    EXPECT_FALSE(water.has_vessel());

    // offsets
    double x_off = 0.0;
    double y_off = 0.0;

    // near/far offsets
    // nearR = 123.634986958 farR = 157.883749639

    // add some cells without bisecting vessel
    bool caught = false;
    try
    {
        x_off = -21.62 * 6.0;
        y_off = 30.0 * 2.0;
        RTK_Cell vessel(10, 21.62, 30.0, 14.28, 1.0, 2.0, x_off, y_off, 101);
    }
    catch(const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        x_off = 21.62 * 5.0;
        y_off = 30.0 * 2.0;
        RTK_Cell vessel(10, 21.62, 30.0, 14.28, 35.0, 36.0, x_off, y_off, 101);
    }
    catch(const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        x_off = 21.62 * 5.0;
        y_off = -30.0 * 3.0;
        RTK_Cell vessel(10, 21.62, 30.0, 14.28, 3.0, 35.0, x_off, y_off, 101);
    }
    catch(const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    double R0 = 0.0, R1 = 0.0, xc = 0.0, yc = 0.0;

    x_off = 21.62 * 5.0;
    y_off = -30.0 * 3.0;

    RTK_Cell vessel1(10, 21.62, 30.0, 14.28, 140.0, 158.0, x_off, y_off, 101);
    EXPECT_TRUE(vessel1.has_vessel());
    EXPECT_TRUE(vessel1.vessel_data(R0, R1, xc, yc));
    EXPECT_EQ(140.0, R0);
    EXPECT_EQ(-1, R1);
    EXPECT_SOFTEQ(118.91, xc, 1.0e-12);
    EXPECT_SOFTEQ(-75.0, yc, 1.0e-12);

    x_off = -21.62 * 6.0;
    y_off = 30.0 * 2.0;

    RTK_Cell vessel2(10, 21.62, 30.0, 14.28, 120.0, 157.0, x_off, y_off, 101);
    EXPECT_TRUE(vessel2.has_vessel());
    EXPECT_TRUE(vessel2.vessel_data(R0, R1, xc, yc));
    EXPECT_EQ(-1, R0);
    EXPECT_EQ(157, R1);
    EXPECT_SOFTEQ(-118.91, xc, 1.0e-12);
    EXPECT_SOFTEQ(75.0, yc, 1.0e-12);

    x_off = 21.62 * 5.0;
    y_off = 30.0 * 2.0;

    RTK_Cell vessel3(10, 21.62, 30.0, 14.28, 130.0, 150.0, x_off, y_off, 101);
    EXPECT_TRUE(vessel3.has_vessel());
    EXPECT_TRUE(vessel3.vessel_data(R0, R1, xc, yc));
    EXPECT_EQ(130.0, R0);
    EXPECT_EQ(150.0, R1);
    EXPECT_SOFTEQ(118.91, xc, 1.0e-12);
    EXPECT_SOFTEQ(75.0, yc, 1.0e-12);

    x_off = -21.62 * 6.0;
    y_off = -30.0 * 3.0;

    RTK_Cell vessel4(10, 21.62, 30.0, 14.28, 130.0, 150.0, x_off, y_off, 101);
    EXPECT_TRUE(vessel4.has_vessel());
    EXPECT_TRUE(vessel4.vessel_data(R0, R1, xc, yc));
    EXPECT_EQ(130.0, R0);
    EXPECT_EQ(150.0, R1);
    EXPECT_SOFTEQ(-118.91, xc, 1.0e-12);
    EXPECT_SOFTEQ(-75.0, yc, 1.0e-12);
}

//---------------------------------------------------------------------------//

TEST(Vessel, Initialize)
{
    // offsets
    double x_off = 21.62 * 5.0;
    double y_off = 30.0  * 2.0;

    RTK_Cell vessel(10, 21.62, 30.0, 14.28, 140.0, 158.0, x_off, y_off, 101);

    // near/far for these offsets
    // nearR = 123.634986958 farR = 157.883749639

    // check simple initialization
    Geo_State state;
    vessel.initialize(Vector(0.0, 0.53, 0.0), state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_EQ(Geo_State::NONE, state.face);
}

//---------------------------------------------------------------------------//
// See support/tstPin_Cell_Vessel.py for reference solutions.

TEST(Vessel, Track_LoR)
{
    // offsets
    double x_off = 21.62 * 5.0;
    double y_off = 30.0  * 2.0;

    // near/far for these offsets
    // nearR = 123.634986958 farR = 157.883749639

    RTK_Cell water(10, 21.62, 30.0, 14.28, 140.0, 158.0, x_off, y_off, 101);

    Geo_State state;
    Vector r, omega;
    double eps = 1.0e-10;

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector(  0.942738909161,   0.021633902593,  -0.332829270666);
    water.initialize(r, state);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.2018173738e+01, eps);
    EXPECT_EQ(Geo_State::PLUS_X, state.exiting_face);

    r     = Vector( -2.52,  -4.20,  11.80);
    omega = Vector(  0.130844809690,  -0.929686805689,  -0.344328446552);
    water.initialize(r, state);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.1616815398e+01, eps);
    EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector( -0.638224253549,  -0.728778909543,  -0.248094947929);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::NONE, state.exiting_face);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.4437568710e-02, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R0_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.next_region);

    r     = Vector( -2.52,  -4.20,  11.80);
    omega = Vector(  0.077270894187,   0.976053400258,  -0.203344458387);
    water.initialize(r, state);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 6.4107157490e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R0_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::VESSEL, state.next_region);
}

//---------------------------------------------------------------------------//
// See support/tstPin_Cell_Vessel.py for reference solutions.

TEST(Vessel, Track_HiR)
{
    // offsets
    double x_off = 21.62 * 5.0;
    double y_off = 30.0  * 2.0;

    // near/far for these offsets
    // nearR = 123.634986958 farR = 157.883749639

    RTK_Cell water(10, 21.62, 30.0, 14.28, 122.0, 140.0, x_off, y_off, 101);

    Geo_State state;
    Vector r, omega;
    double eps = 1.0e-10;

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector(  0.942738909161,   0.021633902593,  -0.332829270666);
    water.initialize(r, state);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.2018173738e+01, eps);
    EXPECT_EQ(Geo_State::PLUS_X, state.exiting_face);

    r     = Vector( -2.52,  -4.20,  11.80);
    omega = Vector(  0.130844809690,  -0.929686805689,  -0.344328446552);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.1616815398e+01, eps);
    EXPECT_EQ(Geo_State::MINUS_Y, state.exiting_face);

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector( -0.638224253549,  -0.728778909543,  -0.248094947929);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(Geo_State::NONE, state.exiting_face);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.4437568710e-02, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R1_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::VESSEL, state.next_region);

    r     = Vector( -2.52,  -4.20,  11.80);
    omega = Vector(  0.077270894187,   0.976053400258,  -0.203344458387);
    water.initialize(r, state);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 6.4107157490e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R1_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.next_region);
}

//---------------------------------------------------------------------------//
// See support/tstPin_Cell_Vessel.py for reference solutions.

TEST(Vessel, Track_LoHiR)
{
    // offsets
    double x_off = 21.62 * 5.0;
    double y_off = 30.0  * 2.0;

    // near/far for these offsets
    // nearR = 123.634986958 farR = 157.883749639

    RTK_Cell water(10, 21.62, 30.0, 14.28, 135.0, 140.0, x_off, y_off, 101);

    Geo_State state;
    Vector r, omega;
    double eps = 1.0e-10;

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector(  0.942738909161,   0.021633902593,  -0.332829270666);
    water.initialize(r, state);
    water.distance_to_boundary(r, omega, state);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.2018173738e+01, eps);
    EXPECT_EQ(Geo_State::PLUS_X, state.exiting_face);

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector( -0.638224253549,  -0.728778909543,  -0.248094947929);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(Geo_State::NONE, state.exiting_face);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.4437568710e-02, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R1_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::VESSEL, state.next_region);

    r     = Vector( -2.52,  -4.20,  11.80);
    omega = Vector(  0.130844809690,  -0.929686805689,  -0.344328446552);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 3.4045225628e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R0_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.next_region);

    r     = Vector( -2.52,  -4.20,  11.80);
    omega = Vector(  0.077270894187,   0.976053400258,  -0.203344458387);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region,  6.4107157490e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R1_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.next_region);
}

//---------------------------------------------------------------------------//
// See support/tstPin_Cell_Vessel.py for reference solutions.

TEST(Vessel, Track_Hi2Lo)
{
    // offsets
    double x_off = 21.62 * 5.0;
    double y_off = 30.0  * 2.0;

    // near/far for these offsets
    // nearR = 123.634986958 farR = 157.883749639

    RTK_Cell water(10, 21.62, 30.0, 14.28, 135.0, 140.0, x_off, y_off, 101);

    Geo_State state;
    Vector r, omega;
    double eps = 1.0e-10;

    r     = Vector( -0.52,  -0.25,  11.80);
    omega = Vector( -0.638224253549,  -0.728778909543,  -0.248094947929);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::NONE, state.exiting_face);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.4437568710e-02, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R1_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::VESSEL, state.next_region);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(10, water.matid(state.region));

    // move the ray
    r[0] += omega[0] * state.dist_to_next_region;
    r[1] += omega[1] * state.dist_to_next_region;
    r[2] += omega[2] * state.dist_to_next_region;

    // cross the surface
    water.cross_surface(state);

    EXPECT_EQ(Geo_State::R1_VESSEL, state.face);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_EQ(0, state.segment);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

    // next segment (to lower next shell)
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 5.3925058590e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R0_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.next_region);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_EQ(101, water.matid(state.region));

    // move the ray
    r[0] += omega[0] * state.dist_to_next_region;
    r[1] += omega[1] * state.dist_to_next_region;
    r[2] += omega[2] * state.dist_to_next_region;

    // cross the surface
    water.cross_surface(state);

    EXPECT_EQ(Geo_State::R0_VESSEL, state.face);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(0, state.segment);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

    // next segment (to low x-face)
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.0715916120e+01, eps);
    EXPECT_EQ(Geo_State::MINUS_X, state.exiting_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(10, water.matid(state.region));

    // move the ray
    r[0] += omega[0] * state.dist_to_next_region;
    r[1] += omega[1] * state.dist_to_next_region;
    r[2] += omega[2] * state.dist_to_next_region;

    EXPECT_SOFTEQ(-21.62*0.5, r[0], eps);
}

//---------------------------------------------------------------------------//
// See support/tstPin_Cell_Vessel.py for reference solutions.

TEST(Vessel, Track_Lo2Hi)
{
    // offsets
    double x_off = 21.62 * 5.0;
    double y_off = 30.0  * 2.0;

    // near/far for these offsets
    // nearR = 123.634986958 farR = 157.883749639

    RTK_Cell water(10, 21.62, 30.0, 14.28, 135.0, 140.0, x_off, y_off, 101);

    Geo_State state;
    Vector r, omega;
    double eps = 1.0e-10;

    r     = Vector(-10.00,   0.25,  11.80);
    omega = Vector(  0.723129219960,   0.666634749651,  -0.180782304990);
    water.initialize(r, state);
    EXPECT_EQ(Geo_State::NONE, state.exiting_face);
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region,  2.6917439842e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R0_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::VESSEL, state.next_region);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(10, water.matid(state.region));

    // move the ray
    r[0] += omega[0] * state.dist_to_next_region;
    r[1] += omega[1] * state.dist_to_next_region;
    r[2] += omega[2] * state.dist_to_next_region;

    // cross the surface
    water.cross_surface(state);

    EXPECT_EQ(Geo_State::R0_VESSEL, state.face);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_EQ(0, state.segment);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

    // next segment (to lower next shell)
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 5.1303903122e+00, eps);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);
    EXPECT_EQ(Geo_State::R1_VESSEL, state.next_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.next_region);
    EXPECT_EQ(Geo_State::VESSEL, state.region);
    EXPECT_EQ(101, water.matid(state.region));

    // move the ray
    r[0] += omega[0] * state.dist_to_next_region;
    r[1] += omega[1] * state.dist_to_next_region;
    r[2] += omega[2] * state.dist_to_next_region;

    // cross the surface
    water.cross_surface(state);

    EXPECT_EQ(Geo_State::R1_VESSEL, state.face);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(0, state.segment);
    EXPECT_EQ(Geo_State::INTERNAL, state.exiting_face);

    // next segment (to high y-face)
    water.distance_to_boundary(r, omega, state);
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.4303925000e+01, eps);
    EXPECT_EQ(Geo_State::PLUS_Y, state.exiting_face);
    EXPECT_EQ(Geo_State::MODERATOR, state.region);
    EXPECT_EQ(10, water.matid(state.region));

    // move the ray
    r[0] += omega[0] * state.dist_to_next_region;
    r[1] += omega[1] * state.dist_to_next_region;
    r[2] += omega[2] * state.dist_to_next_region;

    EXPECT_SOFTEQ(15.0, r[1], eps);
}

//---------------------------------------------------------------------------//
//                        end of tstRTK_Cell.cc
//---------------------------------------------------------------------------//
