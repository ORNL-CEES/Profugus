//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/geometry/test/tstRTK_Functions.cc
 * \author Seth R Johnson
 * \date   Fri Jan 18 20:12:25 2013
 * \brief  RTK_Functions unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "utils/Definitions.hh"

#include "../RTK_Functions.hh"
#include "../RTK_State.hh"

using def::Space_Vector;
using profugus::RTK_State;
using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
class MoveFromOutside : public ::testing::Test
{
    void SetUp()
    {
        state.escaping_face = RTK_State::NONE;
        lower = Space_Vector(1., 2., 3.);
        upper = Space_Vector(4., 5., 6.);
    }

  protected:
    RTK_State    state;
    Space_Vector lower;
    Space_Vector upper;
};

//---------------------------------------------------------------------------//

TEST_F(MoveFromOutside, north)
{
    state.d_r   = Space_Vector(2., 1., 4.);
    state.d_dir = Space_Vector(0., 1., 0.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ( 2., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ( 2., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ( 4., state.d_r[Z], 1.e-12);
    EXPECT_EQ(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, north_miss)
{
    state.d_r   = Space_Vector(100., 1., 4.);
    state.d_dir = Space_Vector(  0., 1., 0.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(   4., state.d_r[X], 1.e-12);
    //EXPECT_SOFTEQ(   5., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(   4., state.d_r[Z], 1.e-12);
    EXPECT_NE(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, south)
{
    state.d_r   = Space_Vector(2.,  10., 4.);
    state.d_dir = Space_Vector(0.,  -1., 0.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(  2., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(  5., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(  4., state.d_r[Z], 1.e-12);
    EXPECT_EQ(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, south_miss)
{
    state.d_r   = Space_Vector(100.,  10., 4.);
    state.d_dir = Space_Vector(  0.,  -1., 0.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(   4., state.d_r[X], 1.e-12);
    //EXPECT_SOFTEQ(   2., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(   4., state.d_r[Z], 1.e-12);
    EXPECT_NE(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, east)
{
    state.d_r   = Space_Vector(0., 3., 4.);
    state.d_dir = Space_Vector(1., 0., 0.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(  1., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(  3., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(  4., state.d_r[Z], 1.e-12);
    EXPECT_EQ(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, west)
{
    state.d_r   = Space_Vector(100., 3., 4.);
    state.d_dir = Space_Vector(-1., 0., 0.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(  4., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(  3., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(  4., state.d_r[Z], 1.e-12);
    EXPECT_EQ(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, up)
{
    state.d_r   = Space_Vector(1., 4., 0.);
    state.d_dir = Space_Vector(0., 0., 1.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(  1., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(  4., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(  3., state.d_r[Z], 1.e-12);
    EXPECT_EQ(RTK_State::NONE, state.escaping_face);
}

TEST_F(MoveFromOutside, down)
{
    state.d_r   = Space_Vector(1., 3., 100.);
    state.d_dir = Space_Vector(0., 0., -1.);

    move_from_outside(lower, upper, state);
    EXPECT_SOFTEQ(  1., state.d_r[X], 1.e-12);
    EXPECT_SOFTEQ(  3., state.d_r[Y], 1.e-12);
    EXPECT_SOFTEQ(  6., state.d_r[Z], 1.e-12);
    EXPECT_EQ(RTK_State::NONE, state.escaping_face);
}

//---------------------------------------------------------------------------//
//                        end of tstRTK_Functions.cc
//---------------------------------------------------------------------------//
