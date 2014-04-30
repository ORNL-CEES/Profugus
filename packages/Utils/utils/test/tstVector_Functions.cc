//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstVector_Functions.cc
 * \author Thomas M. Evans
 * \date   Mon Aug 22 12:33:36 2011
 * \brief  Test for general geometry functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>

#include "utils/Definitions.hh"
#include "utils/Constants.hh"
#include "../Vector_Functions.hh"

using def::Space_Vector;
using def::I;
using def::J;
using def::K;

using profugus::cartesian_vector_transform;
using profugus::vector_magnitude;
using profugus::vector_normalize;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Vector, Magnitude)
{
    // make space-vectors
    Space_Vector v1(-1.1, 2.3, 0.9), v2(2.8, 1.3, -9.1);

    // check magnitudes
    double v1r = vector_magnitude(v1);
    double v2r = vector_magnitude(v2);

    EXPECT_SOFTEQ(v1r, 2.7037011669191546, 1.e-12);
    EXPECT_SOFTEQ(v2r, 9.6093704268281801, 1.e-12);
}

TEST(Vector, DotProduct)
{
    // Make space vectors
    Space_Vector v1(-1.1, 2.3, 0.9), v2(2.8, 1.3, -9.1);

    double dp = dot_product(v1, v2);

    EXPECT_SOFTEQ(dp, -8.28, 1.0e-12);
}

TEST(Vector, Normalize)
{
    // make space-vectors
    Space_Vector v1(-1.1, 2.3, 0.9), v2(2.8, 1.3, -9.1);

    // normalize the vectors
    vector_normalize(v1);
    vector_normalize(v2);

    EXPECT_SOFTEQ(v1[I], -1.1 / 2.7037011669191546, 1.e-12);
    EXPECT_SOFTEQ(v1[J], 2.3 / 2.7037011669191546, 1.e-12);
    EXPECT_SOFTEQ(v1[K], 0.9 / 2.7037011669191546, 1.e-12);

    EXPECT_SOFTEQ(v2[I], 2.8 / 9.6093704268281801, 1.e-12);
    EXPECT_SOFTEQ(v2[J], 1.3 / 9.6093704268281801, 1.e-12);
    EXPECT_SOFTEQ(v2[K], -9.1 / 9.6093704268281801, 1.e-12);

    EXPECT_TRUE(soft_equiv(vector_magnitude(v1), 1.0));
    EXPECT_TRUE(soft_equiv(vector_magnitude(v2), 1.0));
}

TEST(Vector, Transform)
{
    Space_Vector v1(-1.1, 2.3, 0.9);
    vector_normalize(v1);

    // transform through some directions
    double costheta = cos(2.0/3.0);
    double sintheta = sqrt(1.0 - costheta*costheta);
    double phi      = 2.0 * profugus::constants::pi / 3.0;

    double a  = 1.0 / sqrt(1.0 - v1[K] * v1[K]);
    double Ox = v1[I] * costheta + v1[K] * v1[I] * sintheta*cos(phi) * a
                - v1[J] * sintheta * sin(phi) * a;
    double Oy = v1[J] * costheta + v1[K] * v1[J] * sintheta*cos(phi) * a
                + v1[I] * sintheta * sin(phi) * a;
    double Oz = v1[K] * costheta - sintheta*cos(phi) / a;

    cartesian_vector_transform(costheta, phi, v1);

    EXPECT_SOFTEQ(v1[I], Ox, 1.e-12);
    EXPECT_SOFTEQ(v1[J], Oy, 1.e-12);
    EXPECT_SOFTEQ(v1[K], Oz, 1.e-12);

    // transform degenerate vector along y
    Space_Vector v3(0.0, 0.0, -1.0), v4(0.0, 0.0, 1.0);

    cartesian_vector_transform(costheta, phi, v3);
    cartesian_vector_transform(costheta, phi, v4);

    EXPECT_TRUE(soft_equiv(v3[I], sintheta*cos(phi)));
    EXPECT_TRUE(soft_equiv(v3[J], sintheta*sin(phi)));
    EXPECT_SOFTEQ(v3[K], -costheta, 1.e-12);

    EXPECT_TRUE(soft_equiv(v4[I], sintheta*cos(phi)));
    EXPECT_TRUE(soft_equiv(v4[J], sintheta*sin(phi)));
    EXPECT_SOFTEQ(v4[K], costheta, 1.e-12);
}

//---------------------------------------------------------------------------//
//                 end of tstVector_Functions.cc
//---------------------------------------------------------------------------//
