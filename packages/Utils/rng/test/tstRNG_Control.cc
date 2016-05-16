//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/rng/test/tstRNG_Control.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 15:06:03 2014
 * \brief  RNG_Control test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "../RNG_Control.hh"

#include "gtest/utils_gtest.hh"

int seed = 2452423;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(RNG_Control, controller)
{
    typedef profugus::RNG_Control::RNG_t RNG;

    // make a controller
    profugus::RNG_Control control(seed);

    // checks
    EXPECT_EQ(1000000000, control.get_number());
    EXPECT_EQ(2452423, control.get_seed());
    EXPECT_EQ(0, control.get_num());

    // make some random numbers
    RNG r0  = control.rng();
    EXPECT_EQ(1, control.get_num());
    RNG r1  = control.rng();
    EXPECT_EQ(2, control.get_num());
    RNG r2  = control.rng();
    EXPECT_EQ(3, control.get_num());

    RNG rr2 = control.rng(2);
    EXPECT_EQ(3, control.get_num());

    RNG rr1 = control.rng(1);
    EXPECT_EQ(2, control.get_num());

    control.set_num(0);
    RNG rr0 = control.rng();
    EXPECT_EQ(1, control.get_num());

    for (int i = 0; i < 100; i++)
    {
        double rn0  = r0.ran();
        double rrn0 = rr0.ran();
        double rn1  = r1.ran();
        double rrn1 = rr1.ran();
        double rn2  = r2.ran();
        double rrn2 = rr2.ran();

        EXPECT_SOFTEQ(rrn0, rn0, 1.0e-12);
        EXPECT_SOFTEQ(rrn1, rn1, 1.0e-12);
        EXPECT_SOFTEQ(rrn2, rn2, 1.0e-12);

        EXPECT_NE(rrn1, rn0);
        EXPECT_NE(rrn2, rn1);
        EXPECT_NE(rrn0, rn2);
    }

    std::vector<char> pack = r0.pack();
    EXPECT_EQ(control.get_size(), pack.size());
}

//---------------------------------------------------------------------------//
//                 end of tstRNG_Control.cc
//---------------------------------------------------------------------------//
