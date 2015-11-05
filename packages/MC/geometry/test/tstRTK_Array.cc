//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/test/tstRTK_Array.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 21 13:16:21 2010
 * \brief  Unit-Test of RTK_Array class.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <memory>

#include "utils/Definitions.hh"
#include "utils/Vector_Functions.hh"
#include "geometry/Definitions.hh"
#include "../RTK_Cell.hh"
#include "../RTK_Array.hh"
#include "../RTK_Geometry.hh"

using namespace std;

typedef profugus::RTK_Cell::Space_Vector Vector;
typedef profugus::RTK_Cell::Geo_State_t  State;

using def::X;
using def::Y;
using def::Z;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(RTKLevels, Levels)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Sub_Lattice;
    typedef profugus::RTK_Array<Sub_Lattice>        Lattice;
    typedef profugus::RTK_Array<Lattice>            Core;

    EXPECT_EQ(0, Sub_Lattice::calc_level());
    EXPECT_EQ(1, Lattice::calc_level());
    EXPECT_EQ(2, Core::calc_level());
}

//---------------------------------------------------------------------------//
/*
 *     0     1     2     3
 *  |-----|-----|-----|-----|
 * 0.0   0.1   0.2   0.3   0.4
 */
TEST(Bounds, all)
{
    vector<double> x(5, 0.0);
    x[1] = 0.1;
    x[2] = 0.2;
    x[3] = 0.3;
    x[4] = 0.4;

    vector<double>::const_iterator itr;
    int                            i;

    // check bounds for interior cells
    itr = lower_bound(x.begin(), x.end(), 0.05);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(0, i);

    itr = lower_bound(x.begin(), x.end(), 0.15);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(1, i);

    itr = lower_bound(x.begin(), x.end(), 0.25);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(2, i);

    itr = lower_bound(x.begin(), x.end(), 0.35);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(3, i);

    // check bounds for faces
    itr = lower_bound(x.begin(), x.end(), 0.0);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(-1, i);
    EXPECT_EQ(x.begin(), itr);

    itr = lower_bound(x.begin(), x.end(), 0.1);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(0, i);
    EXPECT_EQ(x.begin() + 1, itr);

    itr = lower_bound(x.begin(), x.end(), 0.2);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(1, i);
    EXPECT_EQ(x.begin() + 2, itr);

    itr = lower_bound(x.begin(), x.end(), 0.3);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(2, i);
    EXPECT_EQ(x.begin() + 3, itr);

    itr = lower_bound(x.begin(), x.end(), 0.4);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(3, i);
    EXPECT_EQ(x.begin() + 4, itr);

    // check out of bounds
    itr = lower_bound(x.begin(), x.end(), -1.1);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(-1, i);
    EXPECT_EQ(x.begin(), itr);

    itr = lower_bound(x.begin(), x.end(), 1.5);
    i   = itr - x.begin() - 1;
    EXPECT_EQ(4, i);
    EXPECT_EQ(x.end(), itr);

    // transform test
    vector<int> y(5, 0);
    y[0] = 1; y[1] = 3; y[2] = 2; y[3] = 5; y[4] = 1;
    vector<int> yoff(5, 0);
    transform(&yoff[0], &yoff[0] + 4, y.begin(), &yoff[0] + 1, plus<int>());
    EXPECT_EQ(0, yoff[0]);
    EXPECT_EQ(1, yoff[1]);
    EXPECT_EQ(4, yoff[2]);
    EXPECT_EQ(6, yoff[3]);
    EXPECT_EQ(11, yoff[4]);
}

//---------------------------------------------------------------------------//

class LatticeTest : public ::testing::Test
{
  protected:
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice   ;
    typedef Lattice::SP_Object                      SP_Object ;
    typedef Lattice::Object_t                       Object_t  ;
    typedef shared_ptr<Lattice>                     SP_Lattice;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // make a 3x2x1 lattice with 10 pin-types (only 3 defined)
        d_lat = make_shared<Lattice>(3, 2, 1, 10);
        Lattice& lat = *d_lat;

        // pins 0, 3 and 5
        SP_Object pin0(make_shared<Object_t>(0, 0.62, 10, 1.26, 14.28));
        SP_Object pin3(make_shared<Object_t>(3, 0.45, 10, 1.26, 14.28));
        SP_Object pin5(make_shared<Object_t>(5, 0.54, 10, 1.26, 14.28));

        CHECK(2 == pin0->num_cells());
        CHECK(2 == pin3->num_cells());
        CHECK(2 == pin5->num_cells());

        // assign pins to the lattice
        lat.id(0, 0, 0) = 3;
        lat.id(1, 0, 0) = 3;
        lat.id(2, 0, 0) = 3;
        lat.id(0, 1, 0) = 5;
        lat.id(1, 1, 0) = 0;
        lat.id(2, 1, 0) = 5;

        // assign pins to ids
        lat.assign_object(pin0, 0);
        lat.assign_object(pin3, 3);
        lat.assign_object(pin5, 5);
    }

  protected:
    SP_Lattice d_lat;
};

//---------------------------------------------------------------------------//

TEST_F(LatticeTest, construction)
{
    Lattice& lat = *d_lat;
    EXPECT_EQ(0, lat.level());

    // complete
    EXPECT_TRUE(!lat.completed());
    lat.complete(1.0, 2.0, 3.0);
    EXPECT_TRUE(lat.completed());

    EXPECT_EQ(10, lat.num_objects());
    EXPECT_EQ(6, lat.size());
    EXPECT_EQ(3, lat.size(X));
    EXPECT_EQ(2, lat.size(Y));
    EXPECT_EQ(1, lat.size(Z));

    EXPECT_SOFTEQ( 3.78, lat.pitch(X), 1.e-12);
    EXPECT_SOFTEQ( 2.52, lat.pitch(Y), 1.e-12);
    EXPECT_SOFTEQ(14.28, lat.height(), 1.e-12);

    Vector lower, upper;
    lat.get_extents(lower, upper);
    EXPECT_SOFTEQ(        1., lower[X], 1.e-12);
    EXPECT_SOFTEQ(        2., lower[Y], 1.e-12);
    EXPECT_SOFTEQ(        3., lower[Z], 1.e-12);
    EXPECT_SOFTEQ( 3.78 + 1., upper[X], 1.e-12);
    EXPECT_SOFTEQ( 2.52 + 2., upper[Y], 1.e-12);
    EXPECT_SOFTEQ(14.28 + 3., upper[Z], 1.e-12);

    EXPECT_EQ(12, lat.num_cells());
}

//---------------------------------------------------------------------------//

TEST_F(LatticeTest, initialize1)
{
    Lattice& lat = *d_lat; lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(1.261, 2.44, 12.1);
    lat.initialize(r, state);
    EXPECT_EQ(1, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(1, state.region);
    EXPECT_EQ(10, lat.matid(state));
}

//---------------------------------------------------------------------------//

TEST_F(LatticeTest, initialize2)
{
    Lattice& lat = *d_lat; lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(1.259, 1.27, 1.1);
    lat.initialize(r, state);
    EXPECT_EQ(0, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(1, state.region);
    EXPECT_EQ(10, lat.matid(state));
}

//---------------------------------------------------------------------------//

TEST_F(LatticeTest, initialize3)
{
    Lattice& lat = *d_lat;
    lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(3.560000,   2.239887,   1.300000);
    lat.initialize(r, state);
    EXPECT_EQ(2, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(5, lat.matid(state));
}

//---------------------------------------------------------------------------//

TEST_F(LatticeTest, initialize4)
{
    Lattice& lat = *d_lat; lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(1.570000,   0.931993,   2.700000);
    lat.initialize(r, state);
    EXPECT_EQ(1, state.level_coord[0][X]);
    EXPECT_EQ(0, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(3, lat.matid(state));
}

//---------------------------------------------------------------------------//

TEST_F(LatticeTest, initialize5)
{
    Lattice& lat = *d_lat; lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(1.300000,   2.044919,   3.800000);
    lat.initialize(r, state);
    EXPECT_EQ(1, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(0, lat.matid(state));
}

//---------------------------------------------------------------------------//

// Inside (1, 1, 0) Fuel
TEST_F(LatticeTest, fuel)
{
    Lattice& lat = *d_lat; lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(  2.31,   1.99,  14.20);
    Vector omega( -0.262232636986,  -0.030141682412,  -0.964533837188);
    lat.initialize(r, state);
    lat.distance_to_boundary(r, omega, state);

    EXPECT_EQ(1, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(0, lat.matid(state));
    EXPECT_SOFTEQ(state.dist_to_next_region, 3.9647739003, 1.e-6);
    EXPECT_EQ(State::INTERNAL, state.exiting_face);
    EXPECT_EQ(1, state.next_region);
}

//---------------------------------------------------------------------------//

// Inside (0, 1, 0) Moderator
TEST_F(LatticeTest, moderator)
{
    Lattice& lat = *d_lat; lat.complete(0.0, 0.0, 0.0);

    State state;
    Vector r(  1.00,   2.39,   6.20);
    Vector omega( -0.002072117237,  -0.103605861834,   0.994616273607);

    lat.initialize(r, state);
    lat.distance_to_boundary(r, omega, state);

    EXPECT_EQ(0, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(1, state.region);
    EXPECT_EQ(10, lat.matid(state));
    EXPECT_SOFTEQ(state.dist_to_next_region, 1.0107633165, 1.e-6);
    EXPECT_EQ(State::INTERNAL, state.exiting_face);
    EXPECT_EQ(0, state.next_region);

    // check internal boundary crossing
    lat.cross_surface(r, state);

    EXPECT_EQ(0, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(0, state.region);
    EXPECT_EQ(5, lat.matid(state));
    EXPECT_EQ(State::INTERNAL, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Simple 2x2 lattices, see autodoc/core.fig

TEST(Core, all)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice;
    typedef profugus::RTK_Array<Lattice>            Core;

    // 3x3 core with 3 objects (2 lattice types, one reflector) + water on top
    Core core(3, 3, 2, 3);
    {
        // 2 pin types

        // local id = 0
        Lattice::SP_Object black(make_shared<Lattice::Object_t>(
                                     1, 0.54, 5, 1.26, 14.28));
        EXPECT_EQ(2, black->num_cells());

        // local id = 0
        Lattice::SP_Object purple(make_shared<Lattice::Object_t>(
                                      2, 0.54, 5, 1.26, 14.28));
        EXPECT_EQ(2, purple->num_cells());

        // local id = 0
        Lattice::SP_Object water(make_shared<Lattice::Object_t>(
                                     5, 2.52, 14.28));
        EXPECT_EQ(1, water->num_cells());

        // 3 lattices
        Core::SP_Object lat1(make_shared<Lattice>(2, 2, 1, 1));
        Core::SP_Object lat2(make_shared<Lattice>(2, 2, 1, 1));
        Core::SP_Object lat0(make_shared<Lattice>(1, 1, 1, 1));

        // assign pins to ids
        lat1->assign_object(black,  0);
        lat2->assign_object(purple, 0);
        lat0->assign_object(water, 0);

        // assign pins in the lattices
        lat1->id(0, 0, 0) = 0;
        lat1->id(1, 0, 0) = 0;
        lat1->id(0, 1, 0) = 0;
        lat1->id(1, 1, 0) = 0;

        lat2->id(0, 0, 0) = 0;
        lat2->id(1, 0, 0) = 0;
        lat2->id(0, 1, 0) = 0;
        lat2->id(1, 1, 0) = 0;

        lat0->id(0, 0, 0) = 0;

        // complete
        EXPECT_TRUE(!lat1->completed());
        EXPECT_TRUE(!lat2->completed());
        EXPECT_TRUE(!lat0->completed());

        lat1->complete(0.0, 0.0, 0.0);
        lat2->complete(0.0, 0.0, 0.0);
        lat0->complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(lat1->completed());
        EXPECT_TRUE(lat2->completed());
        EXPECT_TRUE(lat0->completed());

        EXPECT_EQ(2.52, lat1->pitch(X));
        EXPECT_EQ(2.52, lat1->pitch(X));
        EXPECT_EQ(2.52, lat2->pitch(Y));
        EXPECT_EQ(2.52, lat2->pitch(Y));
        EXPECT_EQ(2.52, lat0->pitch(X));
        EXPECT_EQ(2.52, lat0->pitch(Y));

        EXPECT_EQ(1, lat0->num_cells());
        EXPECT_EQ(8, lat1->num_cells());
        EXPECT_EQ(8, lat2->num_cells());

        // assign lattices to core
        core.assign_object(lat1, 1);
        core.assign_object(lat2, 2);
        core.assign_object(lat0, 0);

        // assign ids
        core.id(0, 0, 0) = 1;
        core.id(1, 0, 0) = 2;
        core.id(2, 0, 0) = 0;
        core.id(0, 1, 0) = 2;
        core.id(1, 1, 0) = 1;
        core.id(2, 1, 0) = 0;
        core.id(0, 2, 0) = 0;
        core.id(1, 2, 0) = 0;
        core.id(2, 2, 0) = 0;
        core.id(0, 0, 1) = 0;
        core.id(1, 0, 1) = 0;
        core.id(2, 0, 1) = 0;
        core.id(0, 1, 1) = 0;
        core.id(1, 1, 1) = 0;
        core.id(2, 1, 1) = 0;
        core.id(0, 2, 1) = 0;
        core.id(1, 2, 1) = 0;
        core.id(2, 2, 1) = 0;

        // complete assignment
        core.complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(core.completed());

        EXPECT_EQ(46, core.num_cells());
    }

    // check lattice
    EXPECT_TRUE(soft_equiv(core.pitch(X), 7.56));
    EXPECT_TRUE(soft_equiv(core.pitch(Y), 7.56));
    EXPECT_TRUE(soft_equiv(core.height(), 28.56));

    EXPECT_EQ(18, core.size());
    EXPECT_EQ(3, core.size(X));
    EXPECT_EQ(3, core.size(Y));
    EXPECT_EQ(2, core.size(Z));
    EXPECT_EQ(1, core.level());

    // check tracking

    // initialize some points
    Vector r, omega;
    double eps = 1.0e-6, d;
    State  state;
    {
        r     = Vector(  0.044859500000,   5.638180000000,   7.185140000000);
        omega = Vector(  0.994391000000,   0.099447500000,  -0.036019600000);
        core.initialize(r, state);
        core.distance_to_boundary(r, omega, state);

        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(2, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 2.489101872, eps);
        EXPECT_EQ(State::PLUS_X, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(2, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 2.534214409, eps);
        EXPECT_EQ(State::PLUS_X, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(2, state.level_coord[1][X]);
        EXPECT_EQ(2, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 2.534214409, eps);
        EXPECT_EQ(State::PLUS_X, state.exiting_face);

        core.cross_surface(r, state);

        EXPECT_EQ(3, state.level_coord[1][X]);
        EXPECT_EQ(2, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);
        EXPECT_EQ(State::PLUS_X, state.escaping_face);
    }

    // go from water into a fuel pin
    {
        r     = Vector(  4.202350000000,   2.820900000000,  18.507800000000);
        omega = Vector(  0.098705500000,   0.137387000000,  -0.985587000000);

        core.initialize(r, state);
        core.distance_to_boundary(r, omega, state);

        d = state.dist_to_next_region;

        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(1, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 4.289626385, eps);
        EXPECT_EQ(State::MINUS_Z, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(1, core.matid(state));
        EXPECT_SOFTEQ(d, 1.195583643, eps);
        EXPECT_EQ(State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.next_region);
        EXPECT_EQ(0, state.next_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);
        EXPECT_EQ(0, state.face);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 1.495799820, eps);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(1, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 1.505346029, eps);
        EXPECT_EQ(State::PLUS_X, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(2, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 7.665827372, eps);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        // next step
        core.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(2, state.level_coord[1][X]);
        EXPECT_EQ(2, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(5, core.matid(state));
        EXPECT_SOFTEQ(d, 2.626270607, eps);
        EXPECT_EQ(State::MINUS_Z, state.exiting_face);

        // escape
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];

        core.cross_surface(r, state);

        EXPECT_EQ(2, state.level_coord[1][X]);
        EXPECT_EQ(2, state.level_coord[1][Y]);
        EXPECT_EQ(-1, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(-1, state.level_coord[0][Z]);
        EXPECT_EQ(State::MINUS_Z, state.escaping_face);
    }
}

//---------------------------------------------------------------------------//
// Multilevel core with offsets.
/*

  2x2 assemblies
  3x3 core
  4   axial levels in core

 */

TEST(Core, Offset)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice;
    typedef profugus::RTK_Array<Lattice>            Core;

    Core core(3, 3, 4, 1);
    {
        // 1 pin
        Lattice::SP_Object pin(make_shared<Lattice::Object_t>(
                                   1, 0.5, 5, 1.5, 4.0));
        EXPECT_EQ(2, pin->num_cells());

        // 1 lattice
        Core::SP_Object lat(make_shared<Lattice>(2, 2, 1, 1));
        lat->assign_object(pin,  0);
        lat->id(0, 0, 0) = 0;
        lat->id(1, 0, 0) = 0;
        lat->id(0, 1, 0) = 0;
        lat->id(1, 1, 0) = 0;
        lat->complete(0.0, 0.0, 0.0);

        // assign lattices to core
        core.assign_object(lat, 0);
        for (int k = 0; k < 4; ++k)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int i = 0; i < 3; ++i)
                {
                    core.id(i, j, k) = 0;
                }
            }
        }

        core.complete(1.1, 0.0, -5.0);
    }

    State  state;

    core.initialize(Vector(1.2, 0.2, 4.5), state);
    EXPECT_EQ(0, state.level_coord[1][X]);
    EXPECT_EQ(0, state.level_coord[1][Y]);
    EXPECT_EQ(2, state.level_coord[1][Z]);
    EXPECT_EQ(0, state.level_coord[0][X]);
    EXPECT_EQ(0, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(1, state.region);

    core.initialize(Vector(8.34, 2.3, -3.1), state);
    EXPECT_EQ(2, state.level_coord[1][X]);
    EXPECT_EQ(0, state.level_coord[1][Y]);
    EXPECT_EQ(0, state.level_coord[1][Z]);
    EXPECT_EQ(0, state.level_coord[0][X]);
    EXPECT_EQ(1, state.level_coord[0][Y]);
    EXPECT_EQ(0, state.level_coord[0][Z]);
    EXPECT_EQ(0, state.region);
}

//---------------------------------------------------------------------------//
// Simple 2x2 lattices, see autodoc/core.fig

TEST(Reflecting, all)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice;
    typedef profugus::RTK_Array<Lattice>            Core;

    // 3x3 core with 3 objects (2 lattice types, one reflector) + water on
    // top; low x,y,z reflecting conditions
    Core core(3, 3, 2, 3);
    {
        // 2 pin types

        // local id = 0
        Lattice::SP_Object black(make_shared<Lattice::Object_t>(
                                     1, 0.54, 5, 1.26, 14.28));

        // local id = 0
        Lattice::SP_Object purple(make_shared<Lattice::Object_t>(
                                      2, 0.54, 5, 1.26, 14.28));

        // local id = 0
        Lattice::SP_Object water(make_shared<Lattice::Object_t>(
                                     5, 2.52, 14.28));

        // 3 lattices
        Core::SP_Object lat1(make_shared<Lattice>(2, 2, 1, 1));
        Core::SP_Object lat2(make_shared<Lattice>(2, 2, 1, 1));
        Core::SP_Object lat0(make_shared<Lattice>(1, 1, 1, 1));

        // assign pins to ids
        lat1->assign_object(black,  0);
        lat2->assign_object(purple, 0);
        lat0->assign_object(water, 0);

        // assign pins in the lattices
        lat1->id(0, 0, 0) = 0;
        lat1->id(1, 0, 0) = 0;
        lat1->id(0, 1, 0) = 0;
        lat1->id(1, 1, 0) = 0;

        lat2->id(0, 0, 0) = 0;
        lat2->id(1, 0, 0) = 0;
        lat2->id(0, 1, 0) = 0;
        lat2->id(1, 1, 0) = 0;

        lat0->id(0, 0, 0) = 0;

        // complete
        EXPECT_TRUE(!lat1->completed());
        EXPECT_TRUE(!lat2->completed());
        EXPECT_TRUE(!lat0->completed());

        lat1->complete(0.0, 0.0, 0.0);
        lat2->complete(0.0, 0.0, 0.0);
        lat0->complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(lat1->completed());
        EXPECT_TRUE(lat2->completed());
        EXPECT_TRUE(lat0->completed());

        EXPECT_EQ(2.52, lat1->pitch(X));
        EXPECT_EQ(2.52, lat1->pitch(X));
        EXPECT_EQ(2.52, lat2->pitch(Y));
        EXPECT_EQ(2.52, lat2->pitch(Y));
        EXPECT_EQ(2.52, lat0->pitch(X));
        EXPECT_EQ(2.52, lat0->pitch(Y));

        // assign lattices to core
        core.assign_object(lat1, 1);
        core.assign_object(lat2, 2);
        core.assign_object(lat0, 0);

        // assign ids
        core.id(0, 0, 0) = 1;
        core.id(1, 0, 0) = 2;
        core.id(2, 0, 0) = 0;
        core.id(0, 1, 0) = 2;
        core.id(1, 1, 0) = 1;
        core.id(2, 1, 0) = 0;
        core.id(0, 2, 0) = 0;
        core.id(1, 2, 0) = 0;
        core.id(2, 2, 0) = 0;
        core.id(0, 0, 1) = 0;
        core.id(1, 0, 1) = 0;
        core.id(2, 0, 1) = 0;
        core.id(0, 1, 1) = 0;
        core.id(1, 1, 1) = 0;
        core.id(2, 1, 1) = 0;
        core.id(0, 2, 1) = 0;
        core.id(1, 2, 1) = 0;
        core.id(2, 2, 1) = 0;

        // add reflecting boundary conditions
        Core::Vec_Int reflect(6, 0);
        reflect[0] = 1;
        reflect[2] = 1;
        reflect[4] = 1;
        core.set_reflecting(reflect);

        // complete assignment
        core.complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(core.completed());
        EXPECT_EQ(1, core.level());

        Vector r, omega;
        double eps = 1.0e-6, d;
        State  state;

        // go from water into a fuel pin
        {
            r     = Vector(  4.202350000000,   2.820900000000,  18.507800000000);
            omega = Vector(  0.098705500000,   0.137387000000,  -0.985587000000);

            core.initialize(r, state);
            core.distance_to_boundary(r, omega, state);

            d = state.dist_to_next_region;

            EXPECT_EQ(1, state.level_coord[1][X]);
            EXPECT_EQ(1, state.level_coord[1][Y]);
            EXPECT_EQ(1, state.level_coord[1][Z]);
            EXPECT_EQ(0, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(0, state.region);
            EXPECT_EQ(5, core.matid(state));
            EXPECT_SOFTEQ(d, 4.289626385, eps);
            EXPECT_EQ(State::MINUS_Z, state.exiting_face);

            // process surface
            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            core.cross_surface(r, state);

            // next step
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_EQ(1, state.level_coord[1][X]);
            EXPECT_EQ(1, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(1, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(0, state.region);
            EXPECT_EQ(1, core.matid(state));
            EXPECT_SOFTEQ(d, 1.195583643, eps);
            EXPECT_EQ(State::INTERNAL, state.exiting_face);
            EXPECT_EQ(1, state.next_region);
            EXPECT_EQ(0, state.next_face);

            // process surface
            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            core.cross_surface(r, state);
            EXPECT_EQ(0, state.face);

            // next step
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_EQ(1, state.level_coord[1][X]);
            EXPECT_EQ(1, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(1, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(1, state.region);
            EXPECT_EQ(5, core.matid(state));
            EXPECT_SOFTEQ(d, 1.495799820, eps);
            EXPECT_EQ(State::PLUS_Y, state.exiting_face);

            // process surface
            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            core.cross_surface(r, state);

            // next step
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_EQ(1, state.level_coord[1][X]);
            EXPECT_EQ(1, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(1, state.level_coord[0][X]);
            EXPECT_EQ(1, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(1, state.region);
            EXPECT_EQ(5, core.matid(state));
            EXPECT_SOFTEQ(d, 1.505346029, eps);
            EXPECT_EQ(State::PLUS_X, state.exiting_face);

            // process surface
            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            core.cross_surface(r, state);

            // next step
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_EQ(2, state.level_coord[1][X]);
            EXPECT_EQ(1, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(0, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(0, state.region);
            EXPECT_EQ(5, core.matid(state));
            EXPECT_SOFTEQ(d, 7.665827372, eps);
            EXPECT_EQ(State::PLUS_Y, state.exiting_face);

            // process surface
            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            core.cross_surface(r, state);

            // next step
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_EQ(2, state.level_coord[1][X]);
            EXPECT_EQ(2, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(0, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(0, state.region);
            EXPECT_EQ(5, core.matid(state));
            EXPECT_SOFTEQ(d, 2.626270607, eps);
            EXPECT_EQ(State::MINUS_Z, state.exiting_face);

            // reflecting
            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            core.cross_surface(r, state);

            EXPECT_EQ(2, state.level_coord[1][X]);
            EXPECT_EQ(2, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(0, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);
            EXPECT_EQ(State::MINUS_Z, state.reflecting_face);
            EXPECT_EQ(State::MINUS_Z, state.exiting_face);

            // flip direction and make sure that next-distance-to-boundary
            // clears the reflecting flag
            omega[2] = omega[2] - 2.0 * omega[2];
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_EQ(2, state.level_coord[1][X]);
            EXPECT_EQ(2, state.level_coord[1][Y]);
            EXPECT_EQ(0, state.level_coord[1][Z]);
            EXPECT_EQ(0, state.level_coord[0][X]);
            EXPECT_EQ(0, state.level_coord[0][Y]);
            EXPECT_EQ(0, state.level_coord[0][Z]);

            EXPECT_EQ(0, state.region);
            EXPECT_EQ(5, core.matid(state));
            EXPECT_SOFTEQ(d, 14.488827470, eps);
            EXPECT_EQ(State::NONE, state.reflecting_face);
            EXPECT_EQ(State::PLUS_Z, state.exiting_face);
        }
    }

}

//---------------------------------------------------------------------------//
// See support/lattice.png for diagram of lattice

TEST(Ml, Lattice)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Sub_Lattice;
    typedef profugus::RTK_Array<Sub_Lattice>        Lattice;

    // 1 lattice
    Lattice lat(2, 2, 1, 2);

    {
        // 3 pin types

        // local id = 0
        Sub_Lattice::SP_Object black(make_shared<Sub_Lattice::Object_t>(
                                         1, 0.54, 5, 1.26, 14.28));
        // local id = 1
        Sub_Lattice::SP_Object purple(make_shared<Sub_Lattice::Object_t>(
                                          2, 0.54, 5, 1.26, 14.28 ));
        // local id = 0
        Sub_Lattice::SP_Object green(make_shared<Sub_Lattice::Object_t>(
                                         4, 1.20, 6, 2.52, 14.28 ));

        // 2 sub-lattices
        Lattice::SP_Object lat0(make_shared<Sub_Lattice>(2, 2, 1, 2));
        Lattice::SP_Object lat1(make_shared<Sub_Lattice>(1, 1, 1, 1));

        // assign pins to ids
        lat0->assign_object(black,  0);
        lat0->assign_object(purple, 1);

        lat1->assign_object(green, 0);

        // assign pins to the sub-lattices
        lat0->id(0, 0, 0) = 1;
        lat0->id(1, 0, 0) = 0;
        lat0->id(0, 1, 0) = 0;
        lat0->id(1, 1, 0) = 1;

        lat1->id(0, 0, 0) = 0;

        // complete
        EXPECT_TRUE(!lat0->completed());
        EXPECT_TRUE(!lat1->completed());

        lat0->complete(0.0, 0.0, 0.0);
        lat1->complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(lat0->completed());
        EXPECT_TRUE(lat1->completed());

        EXPECT_EQ(2.52, lat1->pitch(X));
        EXPECT_EQ(2.52, lat1->pitch(Y));

        // assign sub-lattices to lattice
        lat.assign_object(lat0, 0);
        lat.assign_object(lat1, 1);

        // assign ids
        lat.id(0, 0, 0) = 0;
        lat.id(1, 1, 0) = 0;
        lat.id(1, 0, 0) = 1;
        lat.id(0, 1, 0) = 1;

        // complete assignment
        lat.complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(lat.completed());
    }

    // check lattice
    EXPECT_TRUE(soft_equiv(lat.pitch(X), 5.04));
    EXPECT_TRUE(soft_equiv(lat.pitch(Y), 5.04));
    EXPECT_TRUE(soft_equiv(lat.height(), 14.28));

    EXPECT_EQ(4, lat.size());
    EXPECT_EQ(2, lat.size(X));
    EXPECT_EQ(2, lat.size(Y));
    EXPECT_EQ(1, lat.size(Z));
    EXPECT_EQ(1, lat.level());

    // check tracking

    // initialize some points
    Vector r, omega;
    double eps = 1.0e-6, d;
    State  state;
    {
        r     = Vector(  0.37,   0.50,   6.20);
        omega = Vector(  0.444820799653,   0.627277235439,   0.639263424651);
        lat.initialize(r, state);
        lat.distance_to_boundary(r, omega, state);

        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(0, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(2, lat.matid(state));
        EXPECT_SOFTEQ(d, 1.012762262, eps);
        EXPECT_EQ(State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.next_region);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];
        state.region = state.next_region;
        state.face   = 0;

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(0, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(5, lat.matid(state));
        EXPECT_SOFTEQ(d, 0.198823233, eps);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];
        state.region = 1;
        state.level_coord[0][Y] = 1;
        state.face = State::NONE;

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(0, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(1, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(5, lat.matid(state));
        EXPECT_SOFTEQ(d, 0.789220224, eps);
        EXPECT_EQ(State::PLUS_X, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];
        state.region = 1;
        state.level_coord[0][X] = 1;
        state.face = State::NONE;

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(0, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(1, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(5, lat.matid(state));
        EXPECT_SOFTEQ(d, 0.202459958, eps);
        EXPECT_EQ(State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);

        // don't hit next pin, step a bit and change direction
        r[0]  = r[0] + 0.5 * d * omega[0];
        r[1]  = r[1] + 0.5 * d * omega[1];
        r[2]  = r[2] + 0.5 * d * omega[2];
        omega = Vector( -0.022352594210,   0.677082719812,   0.735567367455);

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(0, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(1, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(5, lat.matid(state));
        EXPECT_SOFTEQ(d, 1.035975130, eps);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];
        state.region = 1;
        state.level_coord[0][X] = 0;
        state.level_coord[0][Y] = 0;
        state.level_coord[1][Y] = 1;
        state.face = State::NONE;

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(6, lat.matid(state));
        EXPECT_SOFTEQ(d, 0.088858843, eps);
        EXPECT_EQ(State::INTERNAL, state.exiting_face);
        EXPECT_EQ(0, state.next_region);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];
        state.region = state.next_region;
        state.face = 0;

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(0, state.region);
        EXPECT_EQ(4, lat.matid(state));
        EXPECT_SOFTEQ(d, 3.542210531, eps);
        EXPECT_EQ(State::INTERNAL, state.exiting_face);
        EXPECT_EQ(1, state.next_region);

        // process surface
        r[0] = r[0] + d * omega[0];
        r[1] = r[1] + d * omega[1];
        r[2] = r[2] + d * omega[2];
        state.region = state.next_region;
        state.face = 0;

        lat.distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;

        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);

        EXPECT_EQ(1, state.region);
        EXPECT_EQ(6, lat.matid(state));
        EXPECT_SOFTEQ(d, 0.090780152, eps);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);
    }


    // check boundary crossing functions
    {
        fill(state.level_coord.begin(), state.level_coord.end(), 0);

        r = Vector(3.69, 3.15, 6.4);
        state.exiting_face = State::INTERNAL;

        // internal pin surface
        state.level_coord[0][X] = 1;
        state.level_coord[0][Y] = 0;
        state.level_coord[0][Z] = 0;

        state.level_coord[1][X] = 1;
        state.level_coord[1][Y] = 1;
        state.level_coord[1][Z] = 0;

        state.next_face    = 0;
        state.next_region  = 1;

        lat.cross_surface(r, state);

        EXPECT_EQ(1, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);
        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(0, state.exiting_level[0]);
        EXPECT_EQ(0, state.exiting_level[1]);
        EXPECT_EQ(0, state.face);
        EXPECT_EQ(1, state.region);

        // ----

        r = Vector(2.52, 3.15, 6.4);
        state.exiting_face = State::MINUS_X;

        // internal face in level 0
        state.level_coord[0][X] = 1;
        state.level_coord[0][Y] = 0;
        state.level_coord[0][Z] = 0;

        state.level_coord[1][X] = 1;
        state.level_coord[1][Y] = 1;
        state.level_coord[1][Z] = 0;

        state.next_face    = 3;
        state.next_region  = 3;

        lat.cross_surface(r, state);

        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);
        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(State::MINUS_X, state.exiting_face);
        EXPECT_EQ(0, state.exiting_level[0]);
        EXPECT_EQ(0, state.exiting_level[1]);
        EXPECT_EQ(State::NONE, state.face);
        EXPECT_EQ(1, state.region);

        // ----

        r = Vector(1.33, 2.52, 6.4);
        state.exiting_face = State::PLUS_Y;

        // external face in level 0
        state.level_coord[0][X] = 1;
        state.level_coord[0][Y] = 1;
        state.level_coord[0][Z] = 0;

        // internal face in level 1
        state.level_coord[1][X] = 0;
        state.level_coord[1][Y] = 0;
        state.level_coord[1][Z] = 0;

        state.next_face    = 3;
        state.next_region  = 3;

        lat.cross_surface(r, state);

        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);
        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);
        EXPECT_EQ(State::PLUS_Y, state.exiting_level[0]);
        EXPECT_EQ(0, state.exiting_level[1]);
        EXPECT_EQ(State::NONE, state.face);
        EXPECT_EQ(1, state.region);

        // ----

        // push off surface
        lat.update_state(state);

        EXPECT_EQ(State::NONE, state.face);
        EXPECT_EQ(1, state.region);

        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(0, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);
        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(State::PLUS_Y, state.exiting_face);
        EXPECT_EQ(State::PLUS_Y, state.exiting_level[0]);
        EXPECT_EQ(0, state.exiting_level[1]);

        // ----

        r = Vector(2.54, 0.00, 6.4);
        state.exiting_face = State::MINUS_Y;

        // external face in level 0
        state.level_coord[0][X] = 0;
        state.level_coord[0][Y] = 0;
        state.level_coord[0][Z] = 0;

        // external face in level 1
        state.level_coord[1][X] = 1;
        state.level_coord[1][Y] = 0;
        state.level_coord[1][Z] = 0;

        state.next_face    = 3;
        state.next_region  = 3;

        lat.cross_surface(r, state);

        EXPECT_EQ(0, state.level_coord[0][X]);
        EXPECT_EQ(-1, state.level_coord[0][Y]);
        EXPECT_EQ(0, state.level_coord[0][Z]);
        EXPECT_EQ(1, state.level_coord[1][X]);
        EXPECT_EQ(-1, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
        EXPECT_EQ(State::MINUS_Y, state.exiting_face);
        EXPECT_EQ(State::MINUS_Y, state.exiting_level[0]);
        EXPECT_EQ(State::MINUS_Y, state.exiting_level[1]);
        EXPECT_EQ(State::NONE, state.face);
        EXPECT_EQ(1, state.region);
        EXPECT_EQ(State::MINUS_Y, state.escaping_face);
    }

    // find objects
    {
        r = Vector(1.26, 0.01, 6.4);
        lat.find_object(r, state);
        EXPECT_EQ(0, state.level_coord[1][X]);
        EXPECT_EQ(0, state.level_coord[1][Y]);
        EXPECT_EQ(0, state.level_coord[1][Z]);
    }
}

//---------------------------------------------------------------------------//
// See support/lattice_cells.png for core figure showing particle path.

TEST(Core, Cells)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice;
    typedef profugus::RTK_Array<Lattice>            Core;

    // make the core; 3 object types (lattices and reflector)
    Core core(2, 2, 1, 3);
    {
        // 4 pin types

        // local id = 0
        vector<int>    ids1(2, 0);
        vector<double> r1(2, 0.0);
        ids1[0] = 0; r1[0] = 0.5;  // fuel
        ids1[1] = 5; r1[1] = 0.54; // clad
        Lattice::SP_Object pin1(make_shared<Lattice::Object_t>(
                                    ids1, r1, 10, 1.26, 14.28, 4));
        EXPECT_EQ(12, pin1->num_cells());

        // local id = 1
        vector<int>    ids2(3, 0);
        vector<double> r2(3, 0.0);
        ids2[0] = 0; r2[0] = 0.27; // fuel region 1
        ids2[1] = 0; r2[1] = 0.5;  // fuel region 2
        ids2[2] = 5; r2[2] = 0.54; // clad
        Lattice::SP_Object pin2(make_shared<Lattice::Object_t>(
                                    ids2, r2, 10, 1.26, 14.28, 4));
        EXPECT_EQ(16, pin2->num_cells());

        // local id = 2
        Lattice::SP_Object box(make_shared<Lattice::Object_t>(
                                   10, 0.126, 1.26, 14.28));
        EXPECT_EQ(1, box->num_cells());

        // local id = 0
        Lattice::SP_Object mod(make_shared<Lattice::Object_t>(
                                   10, 1.2317, 2.52, 14.28));
        EXPECT_EQ(1, box->num_cells());

        // 3 lattices
        Core::SP_Object lat2(make_shared<Lattice>(4, 2, 1, 3));
        Core::SP_Object lat1(make_shared<Lattice>(4, 2, 1, 3));
        Core::SP_Object lat0(make_shared<Lattice>(1, 1, 1, 1));

        // assign pins to ids
        lat0->assign_object(mod, 0);
        lat1->assign_object(pin1,  0);
        lat1->assign_object(pin2, 1);
        lat1->assign_object(box, 2);
        lat2->assign_object(pin1,  0);
        lat2->assign_object(pin2, 1);
        lat2->assign_object(box, 2);

        // assign pins in the lattices
        lat0->id(0, 0, 0) = 0;

        lat1->id(0, 0, 0) = 2;
        lat1->id(1, 0, 0) = 1;
        lat1->id(2, 0, 0) = 0;
        lat1->id(3, 0, 0) = 2;

        lat1->id(0, 1, 0) = 2;
        lat1->id(1, 1, 0) = 0;
        lat1->id(2, 1, 0) = 1;
        lat1->id(3, 1, 0) = 2;

        lat2->id(0, 0, 0) = 2;
        lat2->id(1, 0, 0) = 0;
        lat2->id(2, 0, 0) = 1;
        lat2->id(3, 0, 0) = 2;

        lat2->id(0, 1, 0) = 2;
        lat2->id(1, 1, 0) = 1;
        lat2->id(2, 1, 0) = 0;
        lat2->id(3, 1, 0) = 2;

        lat2->complete(0.0, 0.0, 0.0);
        lat1->complete(0.0, 0.0, 0.0);
        lat0->complete(0.0, 0.0, 0.0);

        EXPECT_TRUE(soft_equiv(lat0->pitch(0), 1.2317));
        EXPECT_TRUE(soft_equiv(lat0->pitch(1), 2.52));
        EXPECT_EQ(1, lat0->num_cells());

        EXPECT_TRUE(soft_equiv(lat1->pitch(0), 2.772));
        EXPECT_TRUE(soft_equiv(lat1->pitch(1), 2.52));
        EXPECT_EQ(60, lat1->num_cells());

        EXPECT_TRUE(soft_equiv(lat2->pitch(0), 2.772));
        EXPECT_TRUE(soft_equiv(lat2->pitch(1), 2.52));
        EXPECT_EQ(60, lat1->num_cells());

        // assign lattices to core
        core.assign_object(lat0, 0);
        core.assign_object(lat1, 1);
        core.assign_object(lat2, 2);

        // assign ids
        core.id(0, 0, 0) = 1;
        core.id(1, 0, 0) = 0;
        core.id(0, 1, 0) = 2;
        core.id(1, 1, 0) = 0;

        // complete assignment
        core.complete(0.0, 0.0, 0.0);
        EXPECT_TRUE(core.completed());
        EXPECT_EQ(1, core.level());

        EXPECT_TRUE(soft_equiv(core.pitch(0), 4.0037));
        EXPECT_TRUE(soft_equiv(core.pitch(1), 5.04));
        EXPECT_EQ(122, core.num_cells());
    }

    Vector r, omega;
    double eps = 1.0e-6, d;
    State  state;

    // initial point
    r[0] = 1.0e-12;
    r[1] = 4.8737;
    r[2] = 0.0;

    // direction
    omega[0] = 3.2331;
    omega[1] = -4.8737;
    omega[2] = 0.1;
    profugus::vector_normalize(omega);

    // initialize geometry
    core.initialize(r, state);
    EXPECT_EQ(91, core.cellid(state));
    EXPECT_EQ(10, core.matid(state));

    // track through array
    cout.precision(6);
    cout << endl;
    /*
    cout << "Tracking through multicell core" << endl;
    cout << "===============================" << endl;
    cout << "Starting r   = " << setw(12) << fixed << r[0]
         << setw(12) << r[1] << setw(12) << r[2] << endl;
    cout << "Starting dir = " << setw(12) << fixed << omega[0]
         << setw(12) << omega[1] << setw(12) << omega[2] << endl;
     */

    d = 0.0;

    int s   = 0;
    int c[] = {91, 99, 98, 97, 105, 106, 107, 67, 64, 63, 62, 68, 69, 70, 89,
               50, 49, 48, 47, 55, 51, 52, 53, 54, 19, 29, 60};
    while (state.escaping_face == State::NONE)
    {
        // get distance to next boundary
        core.distance_to_boundary(r, omega, state);

        // add up distance
        d += state.dist_to_next_region;

        // push particle to next boundary
        r[0] += omega[0] * state.dist_to_next_region;
        r[1] += omega[1] * state.dist_to_next_region;
        r[2] += omega[2] * state.dist_to_next_region;

        EXPECT_EQ(c[s++], core.cellid(state))
            << "Failure in region " << state.region
            << " with matid " << core.matid(state)
            << "; r = " << r;

        /*
        cout << "Transporting " << setw(10) << state.dist_to_next_region
             << " in cell " << core.cellid(state) << "(" << state.region
             << "," << state.segment << ")" << " through material "
             << core.matid(state) << " to " << setw(10) << r[0]
             << setw(10) << r[1] << setw(10) << r[2] << endl;
         */

        // cross the surface
        core.cross_surface(r, state);
    }

    EXPECT_EQ(27, s);
    EXPECT_SOFTEQ(r[0], 3.2331, 1.e-12);
    EXPECT_SOFTEQ(r[1], 0.0, 1.e-12);
    EXPECT_SOFTEQ(d, -4.8737 / omega[1], 1.e-12);

    cout << endl;
}

//---------------------------------------------------------------------------//

TEST(Symmetric, all)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice;
    typedef shared_ptr< Lattice >                   SP_Lattice;

    // make a 3x3x1 lattice with 3 pin-types
    SP_Lattice lat( make_shared<Lattice>(3, 3, 1, 3));
    EXPECT_EQ(0, lat->level());

    // pins 0, 1 and 2
    Lattice::SP_Object pin0(make_shared<Lattice::Object_t>(
                                0, 0.6, 10, 1.26, 14.28));
    Lattice::SP_Object pin1(make_shared<Lattice::Object_t>(
                                1, 0.6, 10, 1.26, 14.28));
    Lattice::SP_Object pin2(make_shared<Lattice::Object_t>(
                                2, 0.6, 10, 1.26, 14.28));

    EXPECT_EQ(2, pin0->num_cells());
    EXPECT_EQ(2, pin1->num_cells());
    EXPECT_EQ(2, pin2->num_cells());

    // assign pins to the lattice
    lat->id(0, 0, 0) = 0;
    lat->id(1, 0, 0) = 1;
    lat->id(2, 0, 0) = 0;
    lat->id(0, 1, 0) = 2;
    lat->id(1, 1, 0) = 0;
    lat->id(2, 1, 0) = 1;
    lat->id(0, 2, 0) = 0;
    lat->id(1, 2, 0) = 2;
    lat->id(2, 2, 0) = 0;

    // assign pins to ids
    lat->assign_object(pin0, 0);
    lat->assign_object(pin1, 1);
    lat->assign_object(pin2, 2);

    // complete
    EXPECT_TRUE(!lat->completed());
    lat->complete(0.0, 0.0, 0.0);
    EXPECT_TRUE(lat->completed());

    EXPECT_EQ(9, lat->size());
    EXPECT_EQ(3, lat->size(X));
    EXPECT_EQ(3, lat->size(Y));
    EXPECT_EQ(1, lat->size(Z));

    EXPECT_TRUE(soft_equiv(lat->pitch(X), 3.78));
    EXPECT_TRUE(soft_equiv(lat->pitch(Y), 3.78));
    EXPECT_TRUE(soft_equiv(lat->height(), 14.28));

    EXPECT_EQ(9 * 2, lat->num_cells());
}

//---------------------------------------------------------------------------//
// extents: [0 0 0] to [2 3 5]
// Testing intialization from outside the geometry

class LatticeInitTest : public ::testing::Test
{
  protected:
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice ;
    typedef Lattice::Object_t                       Object_t;
    typedef profugus::RTK_Geometry<Lattice>         Geometry_t;

    typedef Lattice::SP_Object     SP_Object  ;
    typedef shared_ptr<Lattice>    SP_Lattice ;
    typedef shared_ptr<Geometry_t> SP_Geometry;
  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // make a 1x1x1 lattice with 10 pin-types (only 3 defined)
        d_lat = make_shared<Lattice>(1, 1, 1, 1);
        Lattice& lat = *d_lat;

        SP_Object pin(make_shared<Object_t>(0, 2., 3., 5., 1));

        // assign pins to the lattice
        lat.id(0, 0, 0) = 0;

        // assign pins to ids
        lat.assign_object(pin, 0);
        lat.complete(-1.0, -2.0, -3.0);

        d_geom = make_shared<Geometry_t>(d_lat);
    }

  protected:
    SP_Lattice  d_lat;
    SP_Geometry d_geom;
    State       state;
};

//---------------------------------------------------------------------------//
// Traveling up, should hit bottom face

TEST_F(LatticeInitTest, below)
{
    Geometry_t& geo = *d_geom;

    Vector   r(1., -4., 0.);
    Vector dir(0., 1., 0.);
    geo.initialize(r, dir, state);

    EXPECT_SOFTEQ( 1., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(-2., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ( 0., state.d_r[Z], 1.e-12);
    EXPECT_EQ(State::NONE, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Traveling northeast, should hit bottom face

TEST_F(LatticeInitTest, leftbelow)
{
    Geometry_t& geo = *d_geom;

    Vector   r(-2., -4., 0.);
    Vector dir(sqrt(.5), sqrt(.5), 0.);
    geo.initialize(r, dir, state);

    EXPECT_SOFTEQ( 0., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(-2., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ( 0., state.d_r[Z], 1.e-12);
    EXPECT_EQ(State::NONE, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Traveling northeast from southwest corner

TEST_F(LatticeInitTest, left1)
{
    Geometry_t& geo = *d_geom;

    Vector   r(-3., -2.5, 0.);
    Vector dir(sqrt(.5), sqrt(.5), 0.);
    geo.initialize(r, dir, state);

    EXPECT_SOFTEQ(-1. , state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(-0.5, state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ( 0. , state.d_r[Z], 1.e-12);
    EXPECT_EQ(State::NONE, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Traveling northeast from left side, should hit left

TEST_F(LatticeInitTest, left2)
{
    Geometry_t& geo = *d_geom;

    Vector   r(-3., -1.5, 0.);
    Vector dir(sqrt(.5), sqrt(.5), 0.);
    geo.initialize(r, dir, state);

    EXPECT_SOFTEQ(-1., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ( 0.5, state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ( 0. , state.d_r[Z], 1.e-12);
    EXPECT_EQ(State::NONE, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Traveling northeast from left side, should hit left

TEST_F(LatticeInitTest, left_nohit)
{
    Geometry_t& geo = *d_geom;

    Vector   r(-3., 0., 0.);
    Vector dir(sqrt(.5), sqrt(.5), 0.);
    geo.initialize(r, dir, state);

    EXPECT_SOFTEQ(-1., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ( 1., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ( 0., state.d_r[Z], 1.e-12);
    EXPECT_EQ(State::PLUS_Y, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Traveling downsouthwest from top corner, should hit +Z

TEST_F(LatticeInitTest, right)
{
    Geometry_t& geo = *d_geom;

    Vector   r(2., 1.9, 2.9);
    Vector dir(-sqrt(1./3), -sqrt(1./3), -sqrt(1./3));
    geo.initialize(r, dir, state);

    EXPECT_SOFTEQ(1. , state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(0.9, state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(1.9, state.d_r[Z], 1.e-12);
    EXPECT_EQ(State::NONE, state.exiting_face);
}

//---------------------------------------------------------------------------//
// Simple 3x3 lattice into core with gap.

TEST(Core, Gap)
{
    typedef profugus::RTK_Array<profugus::RTK_Cell> Lattice;
    typedef profugus::RTK_Array<Lattice>          Core;

    // 1x1 core with 1 objects (1 lattice with a gap)
    Core core(1, 1, 1, 1);
    {
        // mod  = 5
        // fuel = 1

        // 9 pin types -> 8 pins with gap regions, 1 center pin

        // center pin
        Lattice::SP_Object p11(make_shared<Lattice::Object_t>(
                                   1, 0.54, 5, 1.26, 14.28));

        // Gap vectors
        Lattice::Object_t::Gap_Vector g00(0.1, 0.0, 0.1, 0.0);
        Lattice::Object_t::Gap_Vector g10(0.0, 0.0, 0.1, 0.0);
        Lattice::Object_t::Gap_Vector g20(0.0, 0.1, 0.1, 0.0);
        Lattice::Object_t::Gap_Vector g01(0.1, 0.0, 0.0, 0.0);
        Lattice::Object_t::Gap_Vector g21(0.0, 0.1, 0.0, 0.0);
        Lattice::Object_t::Gap_Vector g02(0.1, 0.0, 0.0, 0.1);
        Lattice::Object_t::Gap_Vector g12(0.0, 0.0, 0.0, 0.1);
        Lattice::Object_t::Gap_Vector g22(0.0, 0.1, 0.0, 0.1);

        // single shell pins
        Lattice::Object_t::Vec_Int ids(1, 1);
        Lattice::Object_t::Vec_Dbl r(1, 0.54);

        // edge pins
        Lattice::SP_Object p00(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g00));
        Lattice::SP_Object p10(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g10));
        Lattice::SP_Object p20(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g20));
        Lattice::SP_Object p01(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g01));
        Lattice::SP_Object p21(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g21));
        Lattice::SP_Object p02(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g02));
        Lattice::SP_Object p12(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g12));
        Lattice::SP_Object p22(make_shared<Lattice::Object_t>(
                                   ids, r, 5, 1.26, 14.28, g22));

        // 1 lattices
        Core::SP_Object lat(make_shared<Lattice>(3, 3, 1, 9));

        // assign pins to ids in the lattice
        lat->assign_object(p00, 0);
        lat->assign_object(p10, 1);
        lat->assign_object(p20, 2);
        lat->assign_object(p01, 3);
        lat->assign_object(p11, 4);
        lat->assign_object(p21, 5);
        lat->assign_object(p02, 6);
        lat->assign_object(p12, 7);
        lat->assign_object(p22, 8);

        // assign pins in the lattices
        lat->id(0, 0, 0) = 0;
        lat->id(1, 0, 0) = 1;
        lat->id(2, 0, 0) = 2;
        lat->id(0, 1, 0) = 3;
        lat->id(1, 1, 0) = 4;
        lat->id(2, 1, 0) = 5;
        lat->id(0, 2, 0) = 6;
        lat->id(1, 2, 0) = 7;
        lat->id(2, 2, 0) = 8;

        lat->complete(0.0, 0.0, 0.0);

        EXPECT_SOFTEQ(3.98, lat->pitch(X), 1.0e-12);
        EXPECT_SOFTEQ(3.98, lat->pitch(Y), 1.0e-12);

        // assign lattices to core
        core.assign_object(lat, 0);
        core.id(0, 0, 0) = 0;

        // complete assignment
        core.complete(0.0, 0.0, 0.0);

        EXPECT_TRUE(core.completed());
        EXPECT_SOFTEQ(3.98, core.pitch(X), 1.0e-12);
        EXPECT_SOFTEQ(3.98, core.pitch(Y), 1.0e-12);
        EXPECT_EQ(1, core.size());
        EXPECT_EQ(1, core.size(X));
        EXPECT_EQ(1, core.size(Y));
        EXPECT_EQ(1, core.size(Z));
        EXPECT_EQ(1, core.level());
    }

    // check tracking

    // initialize some points
    Vector r, omega;
    double eps = 1.0e-6, d;
    State  state;

    // track 1
    {
        // see geometry/test/scripts/tstRTK_CoreGap.py for reference calcs
        double ref[] = {0.237519701, 0.554488796, 0.162434912,
                        0.114479720, 0.902761056, 0.200439852,
                        0.364686076, 0.212683515, 1.075883358,
                        0.419384491};

        r     = Vector(  0.05,   0.60,   3.80);
        omega = Vector(  0.603096313805,   0.796275601821,   0.047116899516);
        core.initialize(r, state);

        int s = 0;

        while (state.escaping_face == State::NONE)
        {
            core.distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            EXPECT_SOFTEQ(ref[s], d, 1.0e-8);

            r[0] = r[0] + d * omega[0];
            r[1] = r[1] + d * omega[1];
            r[2] = r[2] + d * omega[2];

            cout << "Transporting " << fixed << setw(10) << d
                 << " in cell " << setw(4) << core.cellid(state)
                 << " through material " << core.matid(state)
                 << " to " << setw(10) << r[0] << setw(10) << r[1]
                 << setw(10) << r[2] << endl;

            core.cross_surface(r, state);
            ++s;
        }
        EXPECT_EQ(10, s);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstRTK_Array.cc
//---------------------------------------------------------------------------//
