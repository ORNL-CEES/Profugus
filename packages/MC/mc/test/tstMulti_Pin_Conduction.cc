//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstMulti_Pin_Conduction.cc
 * \author Steven Hamilton
 * \date   Thu Aug 09 10:06:43 2018
 * \brief  Tests for class Multi_Pin_Conduction.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Multi_Pin_Conduction.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "../Assembly_Model.hh"

using mc::Multi_Pin_Conduction;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class MultiPinConductionTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    using AM = mc::Assembly_Model;

  protected:
    void SetUp()
    {
        std::vector<double> x_edges = {0.0, 1.26, 2.52, 3.78};
        std::vector<double> y_edges = {0.0, 1.26, 2.52, 3.78};

        d_Nx = x_edges.size()-1;
        d_Ny = y_edges.size()-1;

        std::vector<AM::PIN_TYPE> pin_map =
            {AM::FUEL, AM::FUEL,  AM::FUEL,
             AM::FUEL, AM::GUIDE, AM::FUEL,
             AM::FUEL, AM::FUEL,  AM::FUEL};

        double height = 360.0;

        // Build assembly model
        d_assembly = std::make_shared<AM>(
            pin_map, x_edges, y_edges, height);
        d_assembly->set_fuel_radius(0.4275);
        d_assembly->set_clad_radius(0.475);

        // Solver parameters
        auto params = Teuchos::rcp(new Teuchos::ParameterList("params"));
        params->set("fuel_conductivity", 0.0287);
        params->set("clad_conductivity", 0.215);
        params->set("delta_r_fuel", 0.05);
        params->set("delta_r_clad", 0.02);

        // Conduction solver
        d_Nz = 10;
        std::vector<double> dz(d_Nz, height / static_cast<double>(d_Nz));
        d_solver = std::make_shared<Multi_Pin_Conduction>(
            d_assembly, params, dz);
    }

  protected:
    // >>> DATA
    int d_Nx, d_Ny, d_Nz;
    std::shared_ptr<AM> d_assembly;
    std::shared_ptr<Multi_Pin_Conduction> d_solver;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(MultiPinConductionTest, three_by_three)
{
    // Cosine power shape, 70 kW total per pin
    std::vector<double> pin_power =
        {2835.491846985, 5441.269007846, 7606.226914165, 9154.973551872,
         9962.038679129, 9962.038679129, 9154.973551872, 7606.226914165,
         5441.269007846, 2835.491846985};
    std::vector<double> pin_channel_temp =
        {565, 570, 575, 580, 585, 590, 595, 600, 605, 610};

    std::vector<double> power(d_Nx * d_Ny * d_Nz);
    std::vector<double> channel_temp(d_Nx * d_Ny * d_Nz);
    std::vector<double> fuel_temp(d_Nx * d_Ny * d_Nz);

    for (int iz = 0; iz < d_Nz; ++iz)
    {
        for (int iy = 0; iy < d_Ny; ++iy)
        {
            for (int ix = 0; ix < d_Nx; ++ix)
            {
                int assembly_id = ix + d_Nx * (iy + d_Ny * iz);
                channel_temp[assembly_id] = pin_channel_temp[iz];
                if (d_assembly->pin_type(ix, iy) == AM::FUEL)
                    power[assembly_id] = pin_power[iz];
                else
                    power[assembly_id] = 0.0;
            }
        }
    }

    d_solver->solve(power, channel_temp, fuel_temp);

    std::cout << "Fuel temperature: ";
    for (auto val : fuel_temp)
        std::cout << val << " ";
    std::cout << std::endl;

    // Reference solution
    std::vector<double> ref_temp =
        { 719.12410831453, 865.761997984752, 988.439009542625,
         1077.62165136384, 1126.48994646139, 1131.48994646139,
         1092.62165136384, 1013.43900954263, 900.761997984751,
         764.12410831453};

    for (int iz = 0; iz < d_Nz; ++iz)
    {
        for (int iy = 0; iy < d_Ny; ++iy)
        {
            for (int ix = 0; ix < d_Nx; ++ix)
            {
                int assembly_id = ix + d_Nx * (iy + d_Ny * iz);
                if (d_assembly->pin_type(ix, iy) == AM::FUEL)
                {
                    EXPECT_SOFTEQ(ref_temp[iz],
                                  fuel_temp[assembly_id],
                                  1e-6);
                }
                else
                {
                    EXPECT_SOFT_EQ(pin_channel_temp[iz],
                                   fuel_temp[assembly_id]);
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstMulti_Pin_Conduction.cc
//---------------------------------------------------------------------------//
