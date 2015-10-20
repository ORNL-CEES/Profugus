//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstRNG.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 15:05:58 2014
 * \brief  RNG test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "../RNG.hh"

using profugus::RNG;
using namespace std;

int seed = 493875348;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(RNG, Randomness)
{
    RNG ran1(0);
    RNG ran2(1);

    // simple random check
    double r1 = 0.0, r2 = 0.0;
    for (int i = 0; i < 10000; i++)
    {
        r1 += ran1.ran();
        r2 += ran2.ran();
    }

    cout << endl;
    cout << "Simple randomness check, r = " << fixed << setw(12)
         << setprecision(6) << r1 / 10000.0 << " using 10000 samples (0.5)"
         << endl;
    cout << endl;

    EXPECT_SOFTEQ(0.5, r1/10000, 0.01);
    EXPECT_SOFTEQ(0.5, r2/10000, 0.01);

    RNG ran3(0);

    double r3 = 0.0;
    for (int i = 0; i < 10000; i++)
    {
        r3 += ran3.ran();
    }

    EXPECT_SOFTEQ(r1, r3, 1.0e-12);

    // now check that 1, 3 give equal streams
    double eps = 0.00001;
    for (int i = 0; i < 10; i++)
    {
        EXPECT_SOFTEQ(ran1.ran(), ran3.ran(), eps);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstRNG.cc
//---------------------------------------------------------------------------//
