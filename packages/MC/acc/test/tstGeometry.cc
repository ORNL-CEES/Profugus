//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/test/tstGeometry.cc
 * \author Thomas M. Evans
 * \date   Wed Oct 29 15:45:15 2014
 * \brief  Test for Geometry
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "../Geometry.hh"

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

}

//---------------------------------------------------------------------------//

TEST_F(GeometryTest, ray_tracing)
{

}

//---------------------------------------------------------------------------//
//                 end of tstGeometry.cc
//---------------------------------------------------------------------------//
