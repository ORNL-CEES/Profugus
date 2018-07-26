//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstSingle_Pin_Conduction.cc
 * \author Steven Hamilton
 * \date   Thu Jul 26 10:45:14 2018
 * \brief  Tests for class Single_Pin_Conduction.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Single_Pin_Conduction.hh"

#include "Utils/gtest/utils_gtest.hh"

using mc::Single_Pin_Conduction;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SinglePinConductionTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS

  protected:
    void SetUp()
    {
        auto params = Teuchos::rcp(new Teuchos::ParameterList("params"));
        params->set("tolerance", 1.0e-6);
        params->set("max_iters", 10);
        params->set("verbosity", std::string("high"));
        params->set("fuel_conductivity", 0.0287);
        params->set("clad_conductivity", 0.215);
        params->set("delta_r_fuel", 0.05);
        params->set("delta_r_clad", 0.02);

        std::vector<double> dz(10, 36.0);

        d_solver = std::make_shared<Single_Pin_Conduction>(params, dz);
    }

  protected:
    // >>> DATA
    std::shared_ptr<Single_Pin_Conduction> d_solver;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SinglePinConductionTest, standard)
{
    d_solver->set_fuel_radius(0.4275);
    d_solver->set_clad_radius(0.475);

    // Cosine power shape, 70 kW total
    std::vector<double> power =
        {2835.491846985, 5441.269007846, 7606.226914165, 9154.973551872,
         9962.038679129, 9962.038679129, 9154.973551872, 7606.226914165,
         5441.269007846, 2835.491846985};
    std::vector<double> channel_temp =
        {565, 570, 575, 580, 585, 590, 595, 600, 605, 610};
    std::vector<double> fuel_temp(10);

    d_solver->solve(power, channel_temp, fuel_temp);

    // Reference solution
    std::vector<double> ref_temp =
        { 719.12410831453, 865.761997984752, 988.439009542625,
         1077.62165136384, 1126.48994646139, 1131.48994646139,
         1092.62165136384, 1013.43900954263, 900.761997984751,
         764.12410831453};

    EXPECT_VEC_SOFTEQ(ref_temp, fuel_temp, 1.0e-6);
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstSingle_Pin_Conduction.cc
//---------------------------------------------------------------------------//
