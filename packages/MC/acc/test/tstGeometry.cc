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

#include "rng/RNG_Control.hh"
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
    int num_rays = 1;

    // fill a vector with random numbers
    std::vector<double> rnd(num_rays * 5, 0.0);

    // make random numbers
    profugus::RNG_Control con(235235);
    auto rng = con.rng();
    for (auto &r : rnd)
    {
        r = rng.ran();
    }

    ray_trace(*geo, num_rays, rnd);
}

//---------------------------------------------------------------------------//
//                 end of tstGeometry.cc
//---------------------------------------------------------------------------//
