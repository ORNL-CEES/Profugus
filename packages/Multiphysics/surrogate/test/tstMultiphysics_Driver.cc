//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstMultiphysics_Driver.cc
 * \author Steven Hamilton
 * \date   Fri Aug 10 10:01:23 2018
 * \brief  Tests for class Multiphysics_Driver.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Multiphysics_Driver.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "Utils/utils/Constants.hh"

using mc::Multiphysics_Driver;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class MultiphysicsDriverTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    using Diff   = mc::Two_Group_Diffusion;
    using AM     = mc::Assembly_Model;
    using SP_Driver = std::shared_ptr<Multiphysics_Driver>;


  protected:
    void SetUp()
    {
    }

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(MultiphysicsDriverTest, single_pin)
{
    using profugus::constants::pi;

    // Build assembly model
    auto x_edges =  {0.0, 1.26};
    auto y_edges =  {0.0, 1.26};
    std::vector<double> dz(60, 6.0);

    double height = 360.0;

    std::vector<AM::PIN_TYPE> pins = {AM::FUEL};
    auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);
    assembly->set_fuel_radius(0.4275);
    assembly->set_clad_radius(0.475);
    assembly->set_guide_radius(0.612);

    // Set driver parameters
    auto params = Teuchos::rcp(new Teuchos::ParameterList("params"));
    params->set("max_iters", 20);
    params->set("tolerance", 1e-4);
    params->set("relaxation_factor", 0.35);

    // Subchannel parameters
    auto subchannel_params = Teuchos::sublist(params, "Subchannel");
    subchannel_params->set("tolerance", 1.0e-6);
    subchannel_params->set("max_iters", 10);
    subchannel_params->set("verbosity", std::string("none"));
    double ref_area = 1.26 * 1.26 - pi * 0.475 * 0.475;
    subchannel_params->set("mass_flow_rate", 0.35 / ref_area);
    subchannel_params->set("inlet_temperature", 565.0);
    subchannel_params->set("exit_pressure", 1.52e7);

    // Conduction parameters
    auto conduction_params = Teuchos::sublist(params, "Conduction");
    conduction_params->set("fuel_conductivity", 0.0287);
    conduction_params->set("clad_conductivity", 0.215);
    conduction_params->set("delta_r_fuel", 0.05);
    conduction_params->set("delta_r_clad", 0.02);

    // Boundary conditions for diffusion solver
    std::vector<Diff::BC_TYPE> bcs =
    {Diff::REFLECT, Diff::REFLECT,
        Diff::REFLECT, Diff::REFLECT,
        Diff::VACUUM,  Diff::VACUUM};

    // Build driver
    auto driver = std::make_shared<Multiphysics_Driver>(
            assembly, params, dz, bcs);

    driver->solve();

    std::cout << std::endl;
    std::cout << "Power: ";
    for (auto val : driver->power())
        std::cout << val << " ";
    std::cout << std::endl << std::endl;
    std::cout << "Fuel temp: ";
    for (auto val : driver->fuel_temperature())
        std::cout << val << " ";
    std::cout << std::endl << std::endl;
    std::cout << "Coolant temp: ";
    for (auto val : driver->coolant_temperature())
        std::cout << val << " ";
    std::cout << std::endl << std::endl;
    std::cout << "Coolant density: ";
    for (auto val : driver->coolant_density())
        std::cout << val << " ";
    std::cout << std::endl;
}

TEST_F(MultiphysicsDriverTest, assembly)
{
#if UTILS_DBC > 0
    SKIP_TEST("Performance test, not run in debug.");
#endif

    using profugus::constants::pi;

    // Build assembly model
    std::vector<double> x_edges(18, 0.0);
    std::vector<double> y_edges(18, 0.0);
    for (int i = 0; i < 17; ++i)
    {
        x_edges[i+1] = x_edges[i] + 1.26;
        y_edges[i+1] = y_edges[i] + 1.26;
    }
    std::vector<double> dz(20, 18.0);

    constexpr auto F = AM::FUEL;
    constexpr auto G = AM::GUIDE;

    std::vector<AM::PIN_TYPE> pins =
        {F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
         F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
         F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
         F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F};
    double height = 360.0;

    auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);
    assembly->set_fuel_radius(0.4275);
    assembly->set_clad_radius(0.475);
    assembly->set_guide_radius(0.612);

    // Set driver parameters
    auto params = Teuchos::rcp(new Teuchos::ParameterList("params"));
    params->set("max_iters", 20);
    params->set("tolerance", 1e-4);
    params->set("relaxation_factor", 0.35);

    // Subchannel parameters
    auto subchannel_params = Teuchos::sublist(params, "Subchannel");
    subchannel_params->set("tolerance", 1.0e-6);
    subchannel_params->set("max_iters", 10);
    subchannel_params->set("verbosity", std::string("none"));
    double ref_area = 1.26 * 1.26 - pi * 0.475 * 0.475;
    subchannel_params->set("mass_flow_rate", 0.35 / ref_area);
    subchannel_params->set("inlet_temperature", 565.0);
    subchannel_params->set("exit_pressure", 1.52e7);

    // Conduction parameters
    auto conduction_params = Teuchos::sublist(params, "Conduction");
    conduction_params->set("fuel_conductivity", 0.0287);
    conduction_params->set("clad_conductivity", 0.215);
    conduction_params->set("delta_r_fuel", 0.05);
    conduction_params->set("delta_r_clad", 0.02);

    // Boundary conditions for diffusion solver
    std::vector<Diff::BC_TYPE> bcs =
        {Diff::REFLECT, Diff::REFLECT,
         Diff::REFLECT, Diff::REFLECT,
         Diff::VACUUM,  Diff::VACUUM};

    // Build driver
    auto driver = std::make_shared<Multiphysics_Driver>(
            assembly, params, dz, bcs);

    driver->solve();
}
//---------------------------------------------------------------------------//
// end of MC/mc/test/tstMultiphysics_Driver.cc
//---------------------------------------------------------------------------//
