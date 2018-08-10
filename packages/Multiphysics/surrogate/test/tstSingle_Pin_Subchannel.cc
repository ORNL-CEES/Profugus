//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstSingle_Pin_Subchannel.cc
 * \author Steven Hamilton
 * \date   Wed Jul 25 16:35:25 2018
 * \brief  Tests for class Single_Pin_Subchannel.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Single_Pin_Subchannel.hh"

#include "Utils/gtest/utils_gtest.hh"

using mc::Single_Pin_Subchannel;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SinglePinSubchannelTest : public ::testing::Test
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

        std::vector<double> dz(10, 36.0);

        d_solver = std::make_shared<Single_Pin_Subchannel>(params, dz);
    }

  protected:
    // >>> DATA
    std::shared_ptr<Single_Pin_Subchannel> d_solver;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SinglePinSubchannelTest, standard)
{
    double area = 1.26*1.26 - 3.1415926 * 0.475 * 0.475;
    d_solver->set_channel_area(area);
    d_solver->set_mass_flow_rate(0.35);
    d_solver->set_inlet_temperature(565.0);
    d_solver->set_exit_pressure(1.52e7);

    std::vector<double> power(10, 7000.0);
    std::vector<double> temperature(10);
    std::vector<double> density(10);

    d_solver->solve(power, temperature, density);

    std::cout << "Computed temperature: ";
    for (auto val : temperature)
        std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "Computed density: ";
    for (auto val : density)
        std::cout << val << " ";
    std::cout << std::endl;

    std::vector<double> ref_temp =
        {566.85737, 570.53003, 574.14490, 577.70045, 581.19513,
         584.62739, 587.99571, 591.29853, 594.53431, 597.70152};
    EXPECT_VEC_SOFTEQ(ref_temp, temperature, 1.0e-6);

    std::vector<double> ref_dens =
        {0.74084773, 0.73308312, 0.72511111, 0.71692592, 0.70852217,
         0.69989494, 0.69103975, 0.68195270, 0.67263047, 0.66307039};
    EXPECT_VEC_SOFTEQ(ref_dens, density, 1.0e-5);
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstSingle_Pin_Subchannel.cc
//---------------------------------------------------------------------------//
