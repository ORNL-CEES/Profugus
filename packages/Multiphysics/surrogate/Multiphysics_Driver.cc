//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Multiphysics_Driver.cc
 * \author Steven Hamilton
 * \date   Fri Aug 10 09:00:11 2018
 * \brief  Multiphysics_Driver class definitions.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Multiphysics_Driver.hh"

#include "Utils/harness/Soft_Equivalence.hh"

namespace mc
{
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Multiphysics_Driver::Multiphysics_Driver(SP_Assembly    assembly,
                                         RCP_PL         params,
                                         const Vec_Dbl& dz,
                                         const Vec_BC&  bcs)
    : d_assembly(assembly)
    , d_dz(dz)
{
    REQUIRE(assembly);

    double height = std::accumulate(dz.begin(), dz.end(), 0.0);
    CHECK(profugus::soft_equiv(height, d_assembly->height()));

    // Resize storage for solutions
    d_Nx = d_assembly->num_pins_x();
    d_Ny = d_assembly->num_pins_y();
    d_Nz = d_dz.size();
    d_N  = d_Nx * d_Ny * d_Nz;

    d_power.resize(d_N);
    d_fuel_temperature.resize(d_N);
    d_coolant_temperature.resize(d_N);
    d_coolant_density.resize(d_N);

    d_power_magnitude = params->get("total_power", 70000*d_Nx*d_Ny);
    d_tol       = params->get("tolerance", 1e-4);
    d_max_iters = params->get("max_iters", 25);
    d_damping   = params->get("relaxation_factor", 1.0);

    auto subchannel_params = Teuchos::sublist(params, "Subchannel");
    d_subchannel = std::make_shared<Multi_Pin_Subchannel>(
        assembly, subchannel_params, dz);

    auto conduction_params = Teuchos::sublist(params, "Conduction");
    d_conduction = std::make_shared<Multi_Pin_Conduction>(
        assembly, conduction_params, dz);

    d_diffusion = std::make_shared<Two_Group_Diffusion>(
        assembly, dz, bcs);
}

//---------------------------------------------------------------------------//
// Solve
//---------------------------------------------------------------------------//
void Multiphysics_Driver::solve()
{
    // Set initial guess for power
    std::fill(d_power.begin(),
              d_power.end(),
              d_power_magnitude / static_cast<double>(d_N));

    // Previous iteration power for convergence checking
    Vec_Dbl old_power;

    for (int it = 0; it < d_max_iters; ++it)
    {
        old_power = d_power;

        d_subchannel->solve(d_power, d_coolant_temperature, d_coolant_density);
        d_conduction->solve(d_power, d_coolant_temperature, d_fuel_temperature);
        d_diffusion->solve(d_fuel_temperature, d_coolant_density, d_power);

        // Scale power
        double power_nrm = 0.0;
        for (const auto& val : d_power)
            power_nrm += val;
        for (auto& val : d_power)
            val *= d_power_magnitude / power_nrm;

        // Compute relative difference in power
        double power_diff = 0.0;
        for (int i = 0; i < d_N; ++i)
            power_diff += std::abs(d_power[i] - old_power[i]);
        power_diff /= d_power_magnitude;
        std::cout << "Relative power difference at iteration " << it << " = "
            << power_diff << std::endl;

        // Check for convergence
        if (power_diff < d_tol)
            break;

        // Apply relaxation
        for (int i = 0; i < d_N; ++i)
        {
            d_power[i] = d_damping * d_power[i] +
                         (1.0 - d_damping) * old_power[i];
        }
    }
}

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
// end of MC/mc/Multiphysics_Driver.cc
//---------------------------------------------------------------------------//
