//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/test/tstRTK_Vessel.cc
 * \author Thomas M. Evans
 * \date   Wednesday April 30 10:26:51 2014
 * \brief  RTK_Array tests with a vessel.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <map>
#include <iomanip>
#include <memory>

#include "utils/Definitions.hh"
#include "utils/Constants.hh"
#include "rng/RNG_Control.hh"
#include "../RTK_Cell.hh"
#include "../RTK_Array.hh"

using namespace std;

//---------------------------------------------------------------------------//
// Test fixtures
//---------------------------------------------------------------------------//

class Base : public testing::Test
{
  protected:
    typedef profugus::RTK_Cell             Cell_t;
    typedef profugus::RTK_Array<Cell_t>    Lattice_t;
    typedef profugus::RTK_Array<Lattice_t> Core_t;
    typedef Cell_t::Geo_State_t            Geo_State;
    typedef Cell_t::Space_Vector           Vector;
    typedef def::Vec_Dbl                   Vec_Dbl;
    typedef def::Vec_Int                   Vec_Int;

  protected:
    ~Base() {/*...*/}

  protected:
    int nodes, node;
    int seed;
};

//---------------------------------------------------------------------------//
/*!
  Lattice map:

  o o o o o  4
  o f f f o  3
  o f f f o  2
  o f f f o  1
  o o o o o  0

  0 1 2 3 4

  Dimensions:

 0.0   4.5   6.0   7.5   9.0  13.5
  |     |     |     |     |     |
  | 4.5 | 1.5 | 1.5 | 1.5 | 4.5 |
  |     |     |     |     |     |
     0     1     2     3     4

                 . 0.75  2.25  6.75

*/
class Lattice_Test : public Base
{
  protected:
    typedef Base::Lattice_t     Lattice;
    typedef Lattice::SP_Object  SP_Cell;
    typedef Lattice::Object_t   Cell;
    typedef shared_ptr<Lattice> SP_Lattice;

  protected:
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
        seed  = 1213212 + node;

        // make cells
        SP_Cell side, bottom, fuel, corner;

        fuel   = make_shared<Cell>(1, 0.6, 10, 1.5, 2.5);
        side   = make_shared<Cell>(10, 4.5, 1.5, 2.5);
        bottom = make_shared<Cell>(10, 1.5, 4.5, 2.5);
        corner = make_shared<Cell>(10, 4.5, 2.5);

        lattice = make_shared<Lattice>(5, 5, 4, 5);

        lattice->assign_object(fuel,   1);
        lattice->assign_object(side,   2);
        lattice->assign_object(bottom, 3);
        lattice->assign_object(corner, 4);

        for (int k = 0; k < 4; ++k)
        {
            lattice->id(0, 0, k) = 4;
            lattice->id(1, 0, k) = 3;
            lattice->id(2, 0, k) = 3;
            lattice->id(3, 0, k) = 3;
            lattice->id(4, 0, k) = 4;

            lattice->id(0, 1, k) = 2;
            lattice->id(1, 1, k) = 1;
            lattice->id(2, 1, k) = 1;
            lattice->id(3, 1, k) = 1;
            lattice->id(4, 1, k) = 2;;

            lattice->id(0, 2, k) = 2;
            lattice->id(1, 2, k) = 1;
            lattice->id(2, 2, k) = 1;
            lattice->id(3, 2, k) = 1;
            lattice->id(4, 2, k) = 2;

            lattice->id(0, 3, k) = 2;
            lattice->id(1, 3, k) = 1;
            lattice->id(2, 3, k) = 1;
            lattice->id(3, 3, k) = 1;
            lattice->id(4, 3, k) = 2;

            lattice->id(0, 4, k) = 4;
            lattice->id(1, 4, k) = 3;
            lattice->id(2, 4, k) = 3;
            lattice->id(3, 4, k) = 3;
            lattice->id(4, 4, k) = 4;
        }
    }

  protected:

    SP_Lattice lattice;
};

//---------------------------------------------------------------------------//
/*!
  Core Map

  0 0 0 0  3
  0 F F 0  2
  0 F F 0  1
  0 0 0 0  0

  0 1 2 3

  Lattice F:

  f f  1
  f f  0

  0 1

  Lattice pitch:

  | 0.1 | 1.5 | 1.5 | 0.1 | = 3.2

  Core dimensions:

 0.0   5.0   8.2  11.4  16.4
  |     |     |     |     |
  | 5.0 | 3.2 | 3.2 | 5.0 |
  |     |     |     |     |
     0     1     2     3

              .    3.2   8.2

 */

class Core_Test : public Base
{
  protected:
    typedef Base::Core_t       Core;
    typedef Core::Object_t     Lattice;
    typedef Core::SP_Object    SP_Lattice;
    typedef Lattice::Object_t  Cell;
    typedef Lattice::SP_Object SP_Cell;
    typedef shared_ptr<Core>   SP_Core;
    typedef Cell::Gap_Vector   Gap_Vector;

  protected:
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
        seed  = 1213212 + node;

        // make cells
        SP_Cell side, bottom, f00, f10, f01, f11, corner;

        side   = make_shared<Cell>(10, 5.0, 3.2, 2.5);
        bottom = make_shared<Cell>(10, 3.2, 5.0, 2.5);
        corner = make_shared<Cell>(10, 5.0, 2.5);

        // fuel cells
        Gap_Vector g00(0.1, 0.0, 0.1, 0.0), g10(0.0, 0.1, 0.1, 0.0),
            g01(0.1, 0.0, 0.0, 0.1), g11(0.0, 0.1, 0.0, 0.1);
        Vec_Dbl r(1, 0.6);
        Vec_Int fid(1, 1);

        f00 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g00);
        f10 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g10);
        f01 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g01);
        f11 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g11);

        // lattices
        SP_Lattice lat_f, lat_s, lat_b, lat_c;

        lat_f = make_shared<Lattice>(2, 2, 1, 4);
        lat_f->assign_object(f00, 0);
        lat_f->assign_object(f10, 1);
        lat_f->assign_object(f01, 2);
        lat_f->assign_object(f11, 3);

        lat_f->id(0, 0, 0) = 0;
        lat_f->id(1, 0, 0) = 1;
        lat_f->id(0, 1, 0) = 2;
        lat_f->id(1, 1, 0) = 3;

        lat_s = make_shared<Lattice>(1, 1, 1, 1);
        lat_s->assign_object(side, 0);

        lat_b = make_shared<Lattice>(1, 1, 1, 1);
        lat_b->assign_object(bottom, 0);

        lat_c = make_shared<Lattice>(1, 1, 1, 1);
        lat_c->assign_object(corner, 0);

        lat_f->complete(0.0, 0.0, 0.0);
        lat_s->complete(0.0, 0.0, 0.0);
        lat_b->complete(0.0, 0.0, 0.0);
        lat_c->complete(0.0, 0.0, 0.0);

        // core
        core = make_shared<Core>(4, 4, 3, 4);
        core->assign_object(lat_f, 0);
        core->assign_object(lat_c, 1);
        core->assign_object(lat_s, 2);
        core->assign_object(lat_b, 3);

        for (int k = 0; k < 3; ++k)
        {
            core->id(0, 0, k) = 1;
            core->id(1, 0, k) = 3;
            core->id(2, 0, k) = 3;
            core->id(3, 0, k) = 1;

            core->id(0, 1, k) = 2;
            core->id(1, 1, k) = 0;
            core->id(2, 1, k) = 0;
            core->id(3, 1, k) = 2;

            core->id(0, 2, k) = 2;
            core->id(1, 2, k) = 0;
            core->id(2, 2, k) = 0;
            core->id(3, 2, k) = 2;

            core->id(0, 3, k) = 1;
            core->id(1, 3, k) = 3;
            core->id(2, 3, k) = 3;
            core->id(3, 3, k) = 1;
        }
    }

  protected:

    SP_Core core;
};

//---------------------------------------------------------------------------//
/*!
  Core Map

  0 0 0 0  3
  0 F F 0  2
  0 F F 0  1
  0 0 0 0  0

  0 1 2 3

  Lattice F:

  f f  1
  f f  0

  0 1

  Lattice pitch:

  | 0.1 | 1.5 | 1.5 | 0.1 | = 3.2

  Core dimensions:

 0.0   5.0   8.2  11.4  16.4
  |     |     |     |     |
  | 5.0 | 3.2 | 3.2 | 5.0 |
  |     |     |     |     |
     0     1     2     3

              .    3.2   8.2

 */

class Core_Baffle_Test : public Base
{
  protected:
    typedef Base::Core_t       Core;
    typedef Core::Object_t     Lattice;
    typedef Core::SP_Object    SP_Lattice;
    typedef Lattice::Object_t  Cell;
    typedef Lattice::SP_Object SP_Cell;
    typedef shared_ptr<Core>   SP_Core;
    typedef Cell::Gap_Vector   Gap_Vector;

  protected:
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
        seed  = 1213212 + node;

        // make fuel cells
        SP_Cell f00, f10, f01, f11;

        // fuel cells
        Gap_Vector g00(0.1, 0.0, 0.1, 0.0), g10(0.0, 0.1, 0.1, 0.0),
            g01(0.1, 0.0, 0.0, 0.1), g11(0.0, 0.1, 0.0, 0.1);
        Vec_Dbl r(1, 0.6);
        Vec_Int fid(1, 1);

        f00 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g00);
        f10 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g10);
        f01 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g01);
        f11 = make_shared<Cell>(fid, r, 10, 1.5, 2.5, g11);

        // lattices
        SP_Lattice lat_f, lat_s, lat_b, lat_c;

        lat_f = make_shared<Lattice>(2, 2, 1, 4);
        lat_f->assign_object(f00, 0);
        lat_f->assign_object(f10, 1);
        lat_f->assign_object(f01, 2);
        lat_f->assign_object(f11, 3);

        lat_f->id(0, 0, 0) = 0;
        lat_f->id(1, 0, 0) = 1;
        lat_f->id(0, 1, 0) = 2;
        lat_f->id(1, 1, 0) = 3;

        // side lattice
        {
            SP_Cell c0, c1, c2;
            c0 = make_shared<Cell>(10, 0.1, 1.6, 2.5);
            c1 = make_shared<Cell>(10, 2.0, 1.6, 2.5);
            c2 = make_shared<Cell>(10, 2.9, 1.6, 2.5);
            lat_s = make_shared<Lattice>(3, 2, 1, 3);
            lat_s->assign_object(c0, 0);
            lat_s->assign_object(c1, 1);
            lat_s->assign_object(c2, 2);
            lat_s->id(0,0,0) = 0;
            lat_s->id(1,0,0) = 1;
            lat_s->id(2,0,0) = 2;
            lat_s->id(0,1,0) = 0;
            lat_s->id(1,1,0) = 1;
            lat_s->id(2,1,0) = 2;
        }

        // bottom/top lattice
        {
            SP_Cell c0, c1, c2;
            c0 = make_shared<Cell>(10, 1.6, 0.1, 2.5);
            c1 = make_shared<Cell>(10, 1.6, 2.0, 2.5);
            c2 = make_shared<Cell>(10, 1.6, 2.9, 2.5);
            lat_b = make_shared<Lattice>(2, 3, 1, 3);
            lat_b->assign_object(c0, 0);
            lat_b->assign_object(c1, 1);
            lat_b->assign_object(c2, 2);
            lat_b->id(0,0,0) = 0;
            lat_b->id(1,0,0) = 0;
            lat_b->id(0,1,0) = 1;
            lat_b->id(1,1,0) = 1;
            lat_b->id(0,2,0) = 2;
            lat_b->id(1,2,0) = 2;
        }

        // corner lattices
        {
            SP_Cell c;
            c = make_shared<Cell>(10, 5.0, 5.0, 2.5);
            lat_c = make_shared<Lattice>(1, 1, 1, 1);
            lat_c->assign_object(c, 0);
        }

        lat_f->complete(0.0, 0.0, 0.0);
        lat_s->complete(0.0, 0.0, 0.0);
        lat_b->complete(0.0, 0.0, 0.0);
        lat_c->complete(0.0, 0.0, 0.0);

        // core
        core = make_shared<Core>(4, 4, 3, 4);
        core->assign_object(lat_f, 0);
        core->assign_object(lat_c, 1);
        core->assign_object(lat_s, 2);
        core->assign_object(lat_b, 3);

        for (int k = 0; k < 3; ++k)
        {
            core->id(0, 0, k) = 1;
            core->id(1, 0, k) = 3;
            core->id(2, 0, k) = 3;
            core->id(3, 0, k) = 1;

            core->id(0, 1, k) = 2;
            core->id(1, 1, k) = 0;
            core->id(2, 1, k) = 0;
            core->id(3, 1, k) = 2;

            core->id(0, 2, k) = 2;
            core->id(1, 2, k) = 0;
            core->id(2, 2, k) = 0;
            core->id(3, 2, k) = 2;

            core->id(0, 3, k) = 1;
            core->id(1, 3, k) = 3;
            core->id(2, 3, k) = 3;
            core->id(3, 3, k) = 1;
        }
    }

  protected:

    SP_Core core;
};

//---------------------------------------------------------------------------//

class Track_Check : public Lattice_Test
{
  protected:
    void SetUp()
    {
        Lattice_Test::SetUp();
    }

    void check(Vector    &r,
               Vector    &omega,
               Geo_State &state)
    {
        double d = 0.0;

        // initialize
        lattice->initialize(r, state);

        // Step 1
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 2
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(12, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 3
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 4
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 5
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(1, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 6
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 7
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 8
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(1, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 9
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 10
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 11
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(1, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 12
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 13
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 14
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(12, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_EQ(Geo_State::NONE, state.escaping_face);

        // Step 15
        lattice->distance_to_boundary(r, omega, state);
        d = state.dist_to_next_region;
        EXPECT_EQ(10, lattice->matid(state));

        // move the ray
        r[0] += omega[0] * d;
        r[1] += omega[1] * d;
        r[2] += omega[2] * d;

        // cross the surface
        lattice->cross_surface(r, state);
        EXPECT_NE(Geo_State::NONE, state.escaping_face);
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Lattice_Test, set_vessel)
{
    EXPECT_EQ(5, lattice->num_objects());

    lattice->set_vessel(6.0, 6.5, 101);
    lattice->complete(0.0, 0.0, 0.0);

    EXPECT_TRUE(lattice->has_vessel());

    // adding the vessel makes 16 new cell per axial level (16*4 = 64) cells
    // with a vessel in them + 5 for the original objects
    EXPECT_EQ(69, lattice->num_objects());

    for (int k = 0; k < 4; ++k)
    {
        EXPECT_TRUE(lattice->object(0,  0, k).has_vessel());
        EXPECT_TRUE(lattice->object(1,  0, k).has_vessel());
        EXPECT_TRUE(lattice->object(2,  0, k).has_vessel());
        EXPECT_TRUE(lattice->object(3,  0, k).has_vessel());
        EXPECT_TRUE(lattice->object(4,  0, k).has_vessel());
        EXPECT_TRUE(lattice->object(0,  1, k).has_vessel());
        EXPECT_FALSE(lattice->object(1, 1, k).has_vessel());
        EXPECT_FALSE(lattice->object(2, 1, k).has_vessel());
        EXPECT_FALSE(lattice->object(3, 1, k).has_vessel());
        EXPECT_TRUE(lattice->object(4,  1, k).has_vessel());
        EXPECT_TRUE(lattice->object(0,  2, k).has_vessel());
        EXPECT_FALSE(lattice->object(1, 2, k).has_vessel());
        EXPECT_FALSE(lattice->object(2, 2, k).has_vessel());
        EXPECT_FALSE(lattice->object(3, 2, k).has_vessel());
        EXPECT_TRUE(lattice->object(4,  2, k).has_vessel());
        EXPECT_TRUE(lattice->object(0,  3, k).has_vessel());
        EXPECT_FALSE(lattice->object(1, 3, k).has_vessel());
        EXPECT_FALSE(lattice->object(2, 3, k).has_vessel());
        EXPECT_FALSE(lattice->object(3, 3, k).has_vessel());
        EXPECT_TRUE(lattice->object(4,  3, k).has_vessel());
        EXPECT_TRUE(lattice->object(0,  4, k).has_vessel());
        EXPECT_TRUE(lattice->object(1,  4, k).has_vessel());
        EXPECT_TRUE(lattice->object(2,  4, k).has_vessel());
        EXPECT_TRUE(lattice->object(3,  4, k).has_vessel());
        EXPECT_TRUE(lattice->object(4,  4, k).has_vessel());
    }

    EXPECT_EQ(5,  lattice->id(4, 2, 0));
    EXPECT_EQ(9,  lattice->id(4, 3, 0));
    EXPECT_EQ(13, lattice->id(2, 4, 0));
    EXPECT_EQ(17, lattice->id(3, 4, 0));
    EXPECT_EQ(21, lattice->id(4, 4, 0));
    EXPECT_EQ(25, lattice->id(0, 2, 0));
    EXPECT_EQ(29, lattice->id(0, 3, 0));
    EXPECT_EQ(33, lattice->id(1, 4, 0));
    EXPECT_EQ(37, lattice->id(0, 4, 0));
    EXPECT_EQ(41, lattice->id(4, 1, 0));
    EXPECT_EQ(45, lattice->id(2, 0, 0));
    EXPECT_EQ(49, lattice->id(3, 0, 0));
    EXPECT_EQ(53, lattice->id(4, 0, 0));
    EXPECT_EQ(57, lattice->id(0, 1, 0));
    EXPECT_EQ(61, lattice->id(1, 0, 0));
    EXPECT_EQ(65, lattice->id(0, 0, 0));

    EXPECT_EQ(8,  lattice->id(4, 2, 3));
    EXPECT_EQ(12, lattice->id(4, 3, 3));
    EXPECT_EQ(16, lattice->id(2, 4, 3));
    EXPECT_EQ(20, lattice->id(3, 4, 3));
    EXPECT_EQ(24, lattice->id(4, 4, 3));
    EXPECT_EQ(28, lattice->id(0, 2, 3));
    EXPECT_EQ(32, lattice->id(0, 3, 3));
    EXPECT_EQ(36, lattice->id(1, 4, 3));
    EXPECT_EQ(40, lattice->id(0, 4, 3));
    EXPECT_EQ(44, lattice->id(4, 1, 3));
    EXPECT_EQ(48, lattice->id(2, 0, 3));
    EXPECT_EQ(52, lattice->id(3, 0, 3));
    EXPECT_EQ(56, lattice->id(4, 0, 3));
    EXPECT_EQ(60, lattice->id(0, 1, 3));
    EXPECT_EQ(64, lattice->id(1, 0, 3));
    EXPECT_EQ(68, lattice->id(0, 0, 3));

    for (int k = 0; k < 4; ++k)
    {
        EXPECT_EQ(1, lattice->id(1, 1, k));
        EXPECT_EQ(1, lattice->id(2, 1, k));
        EXPECT_EQ(1, lattice->id(3, 1, k));
        EXPECT_EQ(1, lattice->id(1, 2, k));
        EXPECT_EQ(1, lattice->id(2, 2, k));
        EXPECT_EQ(1, lattice->id(3, 2, k));
        EXPECT_EQ(1, lattice->id(1, 3, k));
        EXPECT_EQ(1, lattice->id(1, 3, k));
        EXPECT_EQ(1, lattice->id(3, 3, k));
    }

    double R0 = 0.0, R1 = 0.0, xc = 0.0, yc = 0.0;
    double eps = 1.0e-12;

    for (int k = 0; k < 4; ++k)
    {
        EXPECT_TRUE(lattice->object(0, 0, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-4.5, xc, eps);
        EXPECT_SOFTEQ(-4.5, yc, eps);

        EXPECT_TRUE(lattice->object(1, 0, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-1.5, xc, eps);
        EXPECT_SOFTEQ(-4.5, yc, eps);

        EXPECT_TRUE(lattice->object(2, 0, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(0.0, xc, eps);
        EXPECT_SOFTEQ(-4.5, yc, eps);

        EXPECT_TRUE(lattice->object(3, 0, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(1.5, xc, eps);
        EXPECT_SOFTEQ(-4.5, yc, eps);

        EXPECT_TRUE(lattice->object(4, 0, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(4.5, xc, eps);
        EXPECT_SOFTEQ(-4.5, yc, eps);

        EXPECT_TRUE(lattice->object(0, 4, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-4.5, xc, eps);
        EXPECT_SOFTEQ( 4.5, yc, eps);

        EXPECT_TRUE(lattice->object(1, 4, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-1.5, xc, eps);
        EXPECT_SOFTEQ( 4.5, yc, eps);

        EXPECT_TRUE(lattice->object(2, 4, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(0.0, xc, eps);
        EXPECT_SOFTEQ( 4.5, yc, eps);

        EXPECT_TRUE(lattice->object(3, 4, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(1.5, xc, eps);
        EXPECT_SOFTEQ( 4.5, yc, eps);

        EXPECT_TRUE(lattice->object(4, 4, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(4.5, xc, eps);
        EXPECT_SOFTEQ( 4.5, yc, eps);

        EXPECT_TRUE(lattice->object(0, 1, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-4.5, xc, eps);
        EXPECT_SOFTEQ(-1.5, yc, eps);

        EXPECT_TRUE(lattice->object(0, 2, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-4.5, xc, eps);
        EXPECT_SOFTEQ( 0.0, yc, eps);

        EXPECT_TRUE(lattice->object(0, 3, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(-4.5, xc, eps);
        EXPECT_SOFTEQ( 1.5, yc, eps);

        EXPECT_TRUE(lattice->object(4, 1, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(4.5, xc, eps);
        EXPECT_SOFTEQ(-1.5, yc, eps);

        EXPECT_TRUE(lattice->object(4, 2, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(4.5, xc, eps);
        EXPECT_SOFTEQ( 0.0, yc, eps);

        EXPECT_TRUE(lattice->object(4, 3, k).vessel_data(R0, R1, xc, yc));
        EXPECT_SOFTEQ(6.0, R0, eps);
        EXPECT_SOFTEQ(6.5, R1, eps);
        EXPECT_SOFTEQ(4.5, xc, eps);
        EXPECT_SOFTEQ( 1.5, yc, eps);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Track_Check, plus_x)
{
    lattice->set_vessel(6.0, 6.5, 12);
    lattice->complete(0.0, 0.0, 0.0);

    Vector    r, omega;
    Geo_State state;
    r[0] = 0.0;
    r[1] = lattice->pitch(1) * 0.5;
    r[2] = 5.0;

    omega[0] = 1.0;
    omega[1] = 0.0;
    omega[2] = 0.0;

    check(r, omega, state);

    EXPECT_EQ(Geo_State::PLUS_X, state.escaping_face);
}

//---------------------------------------------------------------------------//

TEST_F(Track_Check, minus_x)
{
    lattice->set_vessel(6.0, 6.5, 12);
    lattice->complete(0.0, 0.0, 0.0);

    Vector    r, omega;
    Geo_State state;
    r[0] = lattice->pitch(0);
    r[1] = lattice->pitch(1) * 0.5;
    r[2] = 5.0;

    omega[0] = -1.0;
    omega[1] = 0.0;
    omega[2] = 0.0;

    check(r, omega, state);

    EXPECT_EQ(Geo_State::MINUS_X, state.escaping_face);
}

//---------------------------------------------------------------------------//

TEST_F(Track_Check, plus_y)
{
    lattice->set_vessel(6.0, 6.5, 12);
    lattice->complete(0.0, 0.0, 0.0);

    Vector    r, omega;
    Geo_State state;

    r[0] = lattice->pitch(0) * 0.5;
    r[1] = 0.0;
    r[2] = 5.0;

    omega[0] = 0.0;
    omega[1] = 1.0;
    omega[2] = 0.0;

    check(r, omega, state);

    EXPECT_EQ(Geo_State::PLUS_Y, state.escaping_face);
}

//---------------------------------------------------------------------------//

TEST_F(Track_Check, minus_y)
{
    lattice->set_vessel(6.0, 6.5, 12);
    lattice->complete(0.0, 0.0, 0.0);

    Vector    r, omega;
    Geo_State state;

    r[0] = lattice->pitch(0) * 0.5;
    r[1] = lattice->pitch(1);
    r[2] = 5.0;

    omega[0] = 0.0;
    omega[1] = -1.0;
    omega[2] = 0.0;

    check(r, omega, state);

    EXPECT_EQ(Geo_State::MINUS_Y, state.escaping_face);
}

//---------------------------------------------------------------------------//

TEST_F(Lattice_Test, detailed_track)
{
    lattice->set_vessel(6.0, 6.5, 12);
    lattice->complete(0.0, 0.0, 0.0);

    Vector    r, omega;
    Geo_State state;
    double    d;

    r[0] = 0.0;
    r[1] = 0.85;
    r[2] = 5.1;

    omega[0] = 1.0;
    omega[1] = 0.0;
    omega[2] = 0.0;

    // initialize
    lattice->initialize(r, state);

    // Step 1
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(4.0223636606, d, 1.0e-8);
    EXPECT_EQ(10, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 2
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(4.7763633940e-01, d, 1.0e-8);
    EXPECT_EQ(12, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 3
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(1.1591287885, d, 1.0e-8);
    EXPECT_EQ(12, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 4
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(3.4087121146e-01, d, 1.0e-8);
    EXPECT_EQ(10, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 5
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(1.5, d, 1.0e-8);
    EXPECT_EQ(10, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 6
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(3.4087121146e-01, d, 1.0e-8);
    EXPECT_EQ(10, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 7
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(1.1591287885, d, 1.0e-8);
    EXPECT_EQ(12, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 8
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(4.7763633940e-01, d, 1.0e-8);
    EXPECT_EQ(12, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    // Step 9
    lattice->distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;
    EXPECT_SOFTEQ(4.0223636606, d, 1.0e-8);
    EXPECT_EQ(10, lattice->matid(state));

    // move the ray
    r[0] += omega[0] * d;
    r[1] += omega[1] * d;
    r[2] += omega[2] * d;

    // cross the surface
    lattice->cross_surface(r, state);

    EXPECT_EQ(Geo_State::PLUS_X, state.escaping_face);
}

//---------------------------------------------------------------------------//
#ifdef USE_MC

TEST_F(Lattice_Test, check_volumes)
{
    using profugus::constants::pi;

    lattice->set_vessel(6.0, 6.5, 12);
    lattice->complete(0.0, 0.0, 0.0);

    // make a random number generator
    mc::RNG_Control control(seed);
    mc::RNG_Control::RNG rng = control.rng();

    Vector    r, omega;
    Geo_State state;
    double    d;
    double    Vr = lattice->pitch(0) * lattice->pitch(1) * lattice->height();

    // tallies
    vector<double> tally(20, 0.0);
    double tot_path = 0.0;

    // sample points and calculate volumes
    int np = 100000;
    for (int n = 0; n < np; ++n)
    {
        // sample x,y,z
        r[0] = 0.0;
        r[1] = rng.ran() * 13.5;
        r[2] = rng.ran() * 10.0;

        // direction
        omega[0] = 1.0;
        omega[1] = 0.0;
        omega[2] = 0.0;

        // initialize
        lattice->initialize(r, state);

        while (state.escaping_face == Geo_State::NONE)
        {
            // distance-to-boundary
            lattice->distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;

            // tally the pathlength
            tally[lattice->matid(state)] += d;

            // total path
            tot_path += d;

            // move the ray
            r[0] += omega[0] * d;
            r[1] += omega[1] * d;
            r[2] += omega[2] * d;

            // cross the surface
            lattice->cross_surface(r, state);
        }
    }

    profugus::global_sum(&tally[0], 20);
    profugus::global_sum(&tot_path, 1);

    map<int, double> Vref;

    Vref[1]  = 9.0 * (10.0 * pi * 0.6 * 0.6);
    Vref[12] = pi * (6.5*6.5 - 6.0*6.0) * 10.0;
    Vref[10] = Vr - Vref[12] - Vref[1];

    // process tallies
    double Vtot = 0.0;
    for (map<int, double>::const_iterator itr = Vref.begin();
         itr != Vref.end(); ++itr)
    {
        double result = tally[itr->first] / tot_path;
        double V      = result * Vr;
        Vtot         += V;

        double err = fabs(V-Vref[itr->first]) / Vref[itr->first];

        EXPECT_SOFTEQ(Vref[itr->first], V, 5.0e-4);

        if (node == 0)
        {
            cout << setw(4) << itr->first << setw(12) << setprecision(6)
                 << fixed << V << setw(12) << Vref[itr->first]
                 << scientific << setw(16) << err << endl;
        }
    }
}

#endif
//---------------------------------------------------------------------------//

TEST_F(Core_Test, set_vessel)
{
    EXPECT_EQ(4, core->num_objects());

    for (int k = 0; k < 3; ++k)
    {
        EXPECT_EQ(1, core->object(0, 0, k).num_objects());
        EXPECT_EQ(1, core->object(1, 0, k).num_objects());
        EXPECT_EQ(1, core->object(2, 0, k).num_objects());
        EXPECT_EQ(1, core->object(3, 0, k).num_objects());
        EXPECT_EQ(1, core->object(0, 1, k).num_objects());
        EXPECT_EQ(4, core->object(1, 1, k).num_objects());
        EXPECT_EQ(4, core->object(2, 1, k).num_objects());
        EXPECT_EQ(1, core->object(3, 1, k).num_objects());
        EXPECT_EQ(1, core->object(0, 2, k).num_objects());
        EXPECT_EQ(4, core->object(1, 2, k).num_objects());
        EXPECT_EQ(4, core->object(2, 2, k).num_objects());
        EXPECT_EQ(1, core->object(3, 2, k).num_objects());
        EXPECT_EQ(1, core->object(0, 3, k).num_objects());
        EXPECT_EQ(1, core->object(1, 3, k).num_objects());
        EXPECT_EQ(1, core->object(2, 3, k).num_objects());
        EXPECT_EQ(1, core->object(3, 3, k).num_objects());
    }

    core->set_vessel(7.0, 8.0, 101);
    core->complete(0.0, 0.0, 0.0);

    // there are 4 original objects plus 12 * 3 = 36 objects added (12 around
    // the fuel assemblies with 3 axial levels)
    EXPECT_EQ(40, core->num_objects());

    for (int k = 0; k < 3; ++k)
    {
        EXPECT_EQ(2, core->object(0, 0, k).num_objects());
        EXPECT_EQ(2, core->object(1, 0, k).num_objects());
        EXPECT_EQ(2, core->object(2, 0, k).num_objects());
        EXPECT_EQ(2, core->object(3, 0, k).num_objects());
        EXPECT_EQ(2, core->object(0, 1, k).num_objects());
        EXPECT_EQ(4, core->object(1, 1, k).num_objects());
        EXPECT_EQ(4, core->object(2, 1, k).num_objects());
        EXPECT_EQ(2, core->object(3, 1, k).num_objects());
        EXPECT_EQ(2, core->object(0, 2, k).num_objects());
        EXPECT_EQ(4, core->object(1, 2, k).num_objects());
        EXPECT_EQ(4, core->object(2, 2, k).num_objects());
        EXPECT_EQ(2, core->object(3, 2, k).num_objects());
        EXPECT_EQ(2, core->object(0, 3, k).num_objects());
        EXPECT_EQ(2, core->object(1, 3, k).num_objects());
        EXPECT_EQ(2, core->object(2, 3, k).num_objects());
        EXPECT_EQ(2, core->object(3, 3, k).num_objects());

        EXPECT_TRUE(core->object(0, 0, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(1, 0, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(2, 0, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(3, 0, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(0, 1, k).object(0, 0, 0).has_vessel());
        EXPECT_FALSE(core->object(1, 1, k).object(0, 0, 0).has_vessel());
        EXPECT_FALSE(core->object(2, 1, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(3, 1, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(0, 2, k).object(0, 0, 0).has_vessel());
        EXPECT_FALSE(core->object(1, 2, k).object(0, 0, 0).has_vessel());
        EXPECT_FALSE(core->object(2, 2, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(3, 2, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(0, 3, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(1, 3, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(2, 3, k).object(0, 0, 0).has_vessel());
        EXPECT_TRUE(core->object(3, 3, k).object(0, 0, 0).has_vessel());
    }
}

//---------------------------------------------------------------------------//

TEST_F(Core_Baffle_Test, set_vessel)
{
    EXPECT_EQ(4, core->num_objects());

    for (int k = 0; k < 3; ++k)
    {
        EXPECT_EQ(1, core->object(0, 0, k).num_objects());
        EXPECT_EQ(3, core->object(1, 0, k).num_objects());
        EXPECT_EQ(3, core->object(2, 0, k).num_objects());
        EXPECT_EQ(1, core->object(3, 0, k).num_objects());
        EXPECT_EQ(3, core->object(0, 1, k).num_objects());
        EXPECT_EQ(4, core->object(1, 1, k).num_objects());
        EXPECT_EQ(4, core->object(2, 1, k).num_objects());
        EXPECT_EQ(3, core->object(3, 1, k).num_objects());
        EXPECT_EQ(3, core->object(0, 2, k).num_objects());
        EXPECT_EQ(4, core->object(1, 2, k).num_objects());
        EXPECT_EQ(4, core->object(2, 2, k).num_objects());
        EXPECT_EQ(3, core->object(3, 2, k).num_objects());
        EXPECT_EQ(1, core->object(0, 3, k).num_objects());
        EXPECT_EQ(3, core->object(1, 3, k).num_objects());
        EXPECT_EQ(3, core->object(2, 3, k).num_objects());
        EXPECT_EQ(1, core->object(3, 3, k).num_objects());
    }

    core->set_vessel(7.0, 8.0, 101);
    core->complete(0.0, 0.0, 0.0);

    EXPECT_TRUE(core->has_vessel());

    // there are 4 original objects plus 12 * 3 = 36 objects added (12 around
    // the fuel assemblies with 3 axial levels)
    EXPECT_EQ(40, core->num_objects());

    for (int k = 0; k < 3; ++k)
    {
        EXPECT_EQ(5.0, core->object(0,0,k).pitch(0));
        EXPECT_EQ(5.0, core->object(0,0,k).pitch(1));

        EXPECT_EQ(5.0, core->object(0,1,k).pitch(0));
        EXPECT_EQ(3.2, core->object(0,1,k).pitch(1));
        EXPECT_EQ(5.0, core->object(0,2,k).pitch(0));
        EXPECT_EQ(3.2, core->object(0,2,k).pitch(1));
        EXPECT_EQ(5.0, core->object(0,3,k).pitch(0));
        EXPECT_EQ(5.0, core->object(0,3,k).pitch(1));

        EXPECT_EQ(0.1, core->object(0,1,k).object(0,0,0).pitch(0));
        EXPECT_EQ(2.0, core->object(0,1,k).object(1,0,0).pitch(0));
        EXPECT_EQ(2.9, core->object(0,1,k).object(2,0,0).pitch(0));
        EXPECT_EQ(0.1, core->object(0,1,k).object(0,1,0).pitch(0));
        EXPECT_EQ(2.0, core->object(0,1,k).object(1,1,0).pitch(0));
        EXPECT_EQ(2.9, core->object(0,1,k).object(2,1,0).pitch(0));

        EXPECT_EQ(1.6, core->object(0,1,k).object(0,0,0).pitch(1));
        EXPECT_EQ(1.6, core->object(0,1,k).object(1,0,0).pitch(1));
        EXPECT_EQ(1.6, core->object(0,1,k).object(2,0,0).pitch(1));
        EXPECT_EQ(1.6, core->object(0,1,k).object(0,1,0).pitch(1));
        EXPECT_EQ(1.6, core->object(0,1,k).object(1,1,0).pitch(1));
        EXPECT_EQ(1.6, core->object(0,1,k).object(2,1,0).pitch(1));

        EXPECT_EQ(2.5, core->object(0,1,k).object(0,0,0).height());
        EXPECT_EQ(2.5, core->object(0,1,k).object(1,0,0).height());
        EXPECT_EQ(2.5, core->object(0,1,k).object(2,0,0).height());
        EXPECT_EQ(2.5, core->object(0,1,k).object(0,1,0).height());
        EXPECT_EQ(2.5, core->object(0,1,k).object(1,1,0).height());
        EXPECT_EQ(2.5, core->object(0,1,k).object(2,1,0).height());

        EXPECT_FALSE(core->object(0,1,k).object(0,0,0).has_vessel());
        EXPECT_TRUE(core->object(0,1,k).object(1,0,0).has_vessel());
        EXPECT_FALSE(core->object(0,1,k).object(2,0,0).has_vessel());
        EXPECT_FALSE(core->object(0,1,k).object(0,1,0).has_vessel());
        EXPECT_TRUE(core->object(0,1,k).object(1,1,0).has_vessel());
        EXPECT_FALSE(core->object(0,1,k).object(2,1,0).has_vessel());

        EXPECT_EQ(5.0, core->object(3,1,k).pitch(0));
        EXPECT_EQ(3.2, core->object(3,1,k).pitch(1));
        EXPECT_EQ(5.0, core->object(3,2,k).pitch(0));
        EXPECT_EQ(3.2, core->object(3,2,k).pitch(1));
        EXPECT_EQ(5.0, core->object(3,3,k).pitch(0));
        EXPECT_EQ(5.0, core->object(3,3,k).pitch(1));

        EXPECT_EQ(0.1, core->object(3,1,k).object(0,0,0).pitch(0));
        EXPECT_EQ(2.0, core->object(3,1,k).object(1,0,0).pitch(0));
        EXPECT_EQ(2.9, core->object(3,1,k).object(2,0,0).pitch(0));
        EXPECT_EQ(0.1, core->object(3,1,k).object(0,1,0).pitch(0));
        EXPECT_EQ(2.0, core->object(3,1,k).object(1,1,0).pitch(0));
        EXPECT_EQ(2.9, core->object(3,1,k).object(2,1,0).pitch(0));

        EXPECT_EQ(1.6, core->object(3,1,k).object(0,0,0).pitch(1));
        EXPECT_EQ(1.6, core->object(3,1,k).object(1,0,0).pitch(1));
        EXPECT_EQ(1.6, core->object(3,1,k).object(2,0,0).pitch(1));
        EXPECT_EQ(1.6, core->object(3,1,k).object(0,1,0).pitch(1));
        EXPECT_EQ(1.6, core->object(3,1,k).object(1,1,0).pitch(1));
        EXPECT_EQ(1.6, core->object(3,1,k).object(2,1,0).pitch(1));

        EXPECT_FALSE(core->object(3,1,k).object(0,0,0).has_vessel());
        EXPECT_FALSE(core->object(3,1,k).object(1,0,0).has_vessel());
        EXPECT_TRUE(core->object(3,1,k).object(2,0,0).has_vessel());
        EXPECT_FALSE(core->object(3,1,k).object(0,1,0).has_vessel());
        EXPECT_FALSE(core->object(3,1,k).object(1,1,0).has_vessel());
        EXPECT_TRUE(core->object(3,1,k).object(2,1,0).has_vessel());

        double R0, R1, Xc, Yc;
        core->object(3, 1, k).object(2, 0, 0).vessel_data(R0, R1, Xc, Yc);
        EXPECT_EQ(7.0, R0);
        EXPECT_EQ(8.0, R1);
        EXPECT_SOFTEQ(-2.4, Yc, 1.0e-12);
        EXPECT_SOFTEQ(6.75, Xc, 1.0e-12);

        core->object(0, 1, k).object(1, 0, 0).vessel_data(R0, R1, Xc, Yc);
        EXPECT_EQ(7.0, R0);
        EXPECT_EQ(8.0, R1);
        EXPECT_SOFTEQ(-2.4, Yc, 1.0e-12);
        EXPECT_SOFTEQ(-7.1, Xc, 1.0e-12);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Core_Baffle_Test, check_volumes)
{
    using profugus::constants::pi;

    core->set_vessel(7.0, 8.0, 12);
    core->complete(0.0, 0.0, 0.0);

    // make a random number generator
    profugus::RNG_Control control(seed);
    auto rng = control.rng();

    Vector    r, omega;
    Geo_State state;
    double    d;
    double    Vr = core->pitch(0) * core->pitch(1) * core->height();

    // tallies
    vector<double> tally(20, 0.0);
    double tot_path = 0.0;

    // sample points and calculate volumes
    int np = 100000;
    for (int n = 0; n < np; ++n)
    {
        // sample x,y,z
        r[0] = 0.0;
        r[1] = rng.ran() * 16.4;
        r[2] = rng.ran() * 7.5;

        // direction
        omega[0] = 1.0;
        omega[1] = 0.0;
        omega[2] = 0.0;

        // initialize
        core->initialize(r, state);

        while (state.escaping_face == Geo_State::NONE)
        {
            // distance-to-boundary
            core->distance_to_boundary(r, omega, state);
            d = state.dist_to_next_region;
            // tally the pathlength
            tally[core->matid(state)] += d;

            // total path
            tot_path += d;

            // move the ray
            r[0] += omega[0] * d;
            r[1] += omega[1] * d;
            r[2] += omega[2] * d;

            // cross the surface
            core->cross_surface(r, state);
        }
    }

    profugus::global_sum(&tally[0], 20);
    profugus::global_sum(&tot_path, 1);

    map<int, double> Vref;
    map<int, double> err;

    Vref[1]  = 16.0 * (7.5 * pi * 0.6 * 0.6);
    Vref[12] = pi * (8.0*8.0 - 7.0*7.0) * 7.5;
    Vref[10] = Vr - Vref[12] - Vref[1];

    // process tallies
    double Vtot = 0.0;
    for (map<int, double>::const_iterator itr = Vref.begin();
         itr != Vref.end(); ++itr)
    {
        double result = tally[itr->first] / tot_path;
        double V      = result * Vr;
        Vtot         += V;

        err[itr->first] = fabs(V-Vref[itr->first]) / Vref[itr->first];

        EXPECT_SOFTEQ(Vref[itr->first], V, 5.0e-3);

        if (node == 0)
        {
            cout << setw(4) << itr->first << setw(12) << setprecision(6)
                 << fixed << V << setw(12) << Vref[itr->first]
                 << scientific << setw(16) << err[itr->first] << endl;
        }
    }

    EXPECT_LT(err[12], 6.0e-4);
}

//---------------------------------------------------------------------------//
//                 end of tstRTK_Vessel.cc
//---------------------------------------------------------------------------//
