//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstWater_Properties.cc
 * \author Steven Hamilton
 * \date   Wed Jul 25 15:15:23 2018
 * \brief  Tests for class Water_Properties.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Water_Properties.hh"

#include "Utils/gtest/utils_gtest.hh"

using mc::Water_Properties;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(WaterPropertiesTest, density)
{
    // Check against solutions computed using Matlab implementation
    EXPECT_SOFTEQ(Water_Properties::Density(900, 1.5e7), 1.00851075e3, 1e-5);
    EXPECT_SOFTEQ(Water_Properties::Density(200, 1.2e7), 1.00710415e3, 1e-5);
    EXPECT_SOFTEQ(Water_Properties::Temperature(900, 1.5e7), 2.778198e2, 1e-5);
    EXPECT_SOFTEQ(Water_Properties::Temperature(200, 1.2e7), 2.768655e2, 1e-5);
    EXPECT_SOFTEQ(Water_Properties::Enthalpy(600, 1.5e7), 1.498567e6, 1e-5);
    EXPECT_SOFTEQ(Water_Properties::Enthalpy(550, 4.0e6), 1.22228e6,  1e-5);
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstWater_Properties.cc
//---------------------------------------------------------------------------//
