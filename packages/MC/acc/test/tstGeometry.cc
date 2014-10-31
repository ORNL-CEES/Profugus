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
#include <cmath>

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
    typedef acc::Geometry_State            Geo_State;

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

    void normalize(double *r)
    {
        double n = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        r[0] /= n;
        r[1] /= n;
        r[2] /= n;
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

TEST_F(GeometryTest, tracking)
{
#ifndef _OPENACC

    // make a geometry state
    Geo_State state;
    double    dbnd = 0.0;

    double pos[3] = {1.4, 2.5, 4.9};
    double dir[3] = {1.0, 1.0, 1.0};
    normalize(dir);

    // initialize the state
    geo->initialize(pos, dir, state);

    EXPECT_EQ(1, state.ijk[0]);
    EXPECT_EQ(2, state.ijk[1]);
    EXPECT_EQ(4, state.ijk[2]);

    // get distance to boundary
    dbnd = geo->distance_to_boundary(state);
    EXPECT_SOFTEQ(0.173205080757, dbnd, 1.0e-6);

    // step to the boundary
    geo->move_to_surface(state);

    // get distance to boundary
    dbnd = geo->distance_to_boundary(state);
    EXPECT_SOFTEQ(0.692820323028, dbnd, 1.0e-6);

    // step to the boundary
    geo->move_to_surface(state);

    // get distance to boundary
    dbnd = geo->distance_to_boundary(state);
    EXPECT_SOFTEQ(0.173205080757, dbnd, 1.0e-6);

    // step to the boundary
    geo->move_to_surface(state);

    // continue till hitting edge
    while (geo->boundary_state(state) == profugus::geometry::INSIDE)
    {
        dbnd = geo->distance_to_boundary(state);
        geo->move_to_surface(state);
    }

    EXPECT_EQ(profugus::geometry::REFLECT, geo->boundary_state(state));

    EXPECT_SOFTEQ(6.5,  state.pos[0], 1.0e-6);
    EXPECT_SOFTEQ(7.6,  state.pos[1], 1.0e-6);
    EXPECT_SOFTEQ(10.0, state.pos[2], 1.0e-6);

    EXPECT_EQ(6, state.ijk[0]);
    EXPECT_EQ(7, state.ijk[1]);
    EXPECT_EQ(9, state.ijk[2]);

    // reflect the ray
    geo->reflect(state);

    EXPECT_SOFTEQ(0.5773502691896258,  state.dir[0], 1.0e-6);
    EXPECT_SOFTEQ(0.5773502691896258,  state.dir[1], 1.0e-6);
    EXPECT_SOFTEQ(-0.5773502691896258, state.dir[2], 1.0e-6);

    // get distance to boundary
    dbnd = geo->distance_to_boundary(state);
    EXPECT_SOFTEQ(0.692820323028, dbnd, 1.0e-6);

    // step to the boundary
    geo->move_to_surface(state);

#endif
}

//---------------------------------------------------------------------------//

TEST_F(GeometryTest, inf_med)
{
    int num_rays  = 100000;
    int num_steps = 100;

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

    double xl = 0.0;
    for (int i = 0; i < 1000000; ++i)
    {
        xl += -std::log(rng.ran());
    }
    std::cout << "Average log = " << xl / 1000000.0 << std::endl;
}

//---------------------------------------------------------------------------//
//                 end of tstGeometry.cc
//---------------------------------------------------------------------------//
