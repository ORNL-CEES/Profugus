//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstMulti_Pin_Subchannel.cc
 * \author Steven Hamilton
 * \date   Thu Aug 09 13:17:11 2018
 * \brief  Tests for class Multi_Pin_Subchannel.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Multi_Pin_Subchannel.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "Utils/utils/Constants.hh"

#include "../Assembly_Model.hh"

using mc::Multi_Pin_Subchannel;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class MultiPinSubchannelTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    using AM = mc::Assembly_Model;

  protected:
    void SetUp()
    {
        using profugus::constants::pi;

        auto params = Teuchos::rcp(new Teuchos::ParameterList("params"));
        params->set("tolerance", 1.0e-6);
        params->set("max_iters", 10);
        params->set("verbosity", std::string("high"));
        double ref_area = 1.26 * 1.26 - pi * 0.475 * 0.475;
        params->set("mass_flow_rate", 0.35 / ref_area);
        params->set("inlet_temperature", 565.0);
        params->set("exit_pressure", 1.52e7);


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
        d_assembly->set_guide_radius(0.612);

        d_Nz = 10;
        std::vector<double> dz(10, height / static_cast<double>(d_Nz));
        d_solver = std::make_shared<Multi_Pin_Subchannel>(
            d_assembly, params, dz);
    }

  protected:
    // >>> DATA
    int d_Nx, d_Ny, d_Nz;
    std::shared_ptr<AM> d_assembly;
    std::shared_ptr<Multi_Pin_Subchannel> d_solver;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(MultiPinSubchannelTest, three_by_three)
{
    // Cosine power shape, 70 kW total per pin
    std::vector<double> pin_power =
        {2835.491846985, 5441.269007846, 7606.226914165, 9154.973551872,
         9962.038679129, 9962.038679129, 9154.973551872, 7606.226914165,
         5441.269007846, 2835.491846985};

    std::vector<double> power(d_Nx * d_Ny * d_Nz);
    std::vector<double> coolant_temp(d_Nx * d_Ny * d_Nz);
    std::vector<double> coolant_density(d_Nx * d_Ny * d_Nz);

    for (int iz = 0; iz < d_Nz; ++iz)
    {
        for (int iy = 0; iy < d_Ny; ++iy)
        {
            for (int ix = 0; ix < d_Nx; ++ix)
            {
                int assembly_id = ix + d_Nx * (iy + d_Ny * iz);
                if (d_assembly->pin_type(ix, iy) == AM::FUEL)
                    power[assembly_id] = pin_power[iz];
                else
                    power[assembly_id] = 0.0;
            }
        }
    }

    d_solver->solve(power, coolant_temp, coolant_density);

    std::cout << "Coolant temperature: ";
    for (auto val : coolant_temp)
        std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "Coolant density: ";
    for (auto val : coolant_density)
        std::cout << val << " ";
    std::cout << std::endl;
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstMulti_Pin_Subchannel.cc
//---------------------------------------------------------------------------//
