//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/test/tstGeometry.cc
 * \author Thomas M. Evans
 * \date   Wed Oct 29 15:45:15 2014
 * \brief  Test for Geometry
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <memory>
#include <algorithm>

#include "../RNG.hh"
#include "../Geometry.hh"
#include "ParticleTest.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class GeometryTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef std::shared_ptr<acc::Geometry> SP_Geometry;

  protected:
    void SetUp()
    {
        // make geometry
        std::vector<int> bnd = {1, 1, 1, 1, 1, 1};
        std::vector<int> matids(1000, 0);
        for (int i = 1; i < 1000; ++i)
        {
            matids[i] = matids[i-1] + 1;
        }
        geo = std::make_shared<acc::Geometry>(10, 1.0, matids, &bnd[0]);
    }

  protected:
    // >>> DATA

    SP_Geometry geo;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(GeometryTest, construction)
{
    EXPECT_EQ(1000, geo->num_cells());
}

//---------------------------------------------------------------------------//

TEST_F(GeometryTest, ray_tracing)
{
    int num_rays  = 100000;
    int num_steps = 1000;

    // make tallies
    std::vector<double> tallies(geo->num_cells(), 0.0);

    ray_trace(*geo, num_rays, num_steps, tallies);

    // print max and min talliese
    auto maxitr = std::max_element(tallies.begin(), tallies.end());
    auto minitr = std::min_element(tallies.begin(), tallies.end());
    cout << "Maximum pathlength = " << *maxitr << " in "
         << maxitr - tallies.begin() << endl;
    cout << "Minimum pathlength = " << *minitr << " in "
         << minitr - tallies.begin() << endl;
}

//---------------------------------------------------------------------------//

TEST(RNG, randomness)
{
    acc::RNG rng(32534);

    for (int i = 0; i < 10; ++i)
    {
        double total   = 0.0;
        int    samples = 100000;
        for (int i = 0; i < samples; ++i)
        {
            total += rng.ran();
        }

        EXPECT_SOFTEQ(0.5, total / samples, 1.0e-2);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstGeometry.cc
//---------------------------------------------------------------------------//
