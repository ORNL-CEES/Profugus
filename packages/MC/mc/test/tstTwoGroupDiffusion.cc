//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstTwoGroupDiffusion.cc
 * \author Steven Hamilton
 * \date   Tue Aug 07 14:34:50 2018
 * \brief  Tests for class TwoGroupDiffusion.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>

#include "../TwoGroupDiffusion.hh"
#include "../TwoGroupCrossSections.hh"

#include "Utils/gtest/utils_gtest.hh"

using mc::TwoGroupDiffusion;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class TwoGroupDiffusionTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    using Diff = TwoGroupDiffusion;
    using XS   = mc::TwoGroupCrossSections;
    using AM   = mc::Assembly_Model;

  protected:
    void SetUp()
    {
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(TwoGroupDiffusionTest, infhommed)
{
    auto x_edges =  {0.0, 1.26, 2.52, 3.78};
    auto y_edges =  {0.0, 1.26, 2.52, 3.78};
    std::vector<double> dz(2,1.26);
    dz[0] = 1.0;
    dz[1] = 4.0;

    double height = 360.0;

    std::vector<AM::PIN_TYPE> pins =
        {AM::FUEL, AM::FUEL, AM::FUEL,
         AM::FUEL, AM::FUEL, AM::FUEL,
         AM::FUEL, AM::FUEL, AM::FUEL};
    auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);

    std::vector<Diff::BC_TYPE> bcs =
        {Diff::REFLECT, Diff::REFLECT,
         Diff::REFLECT, Diff::REFLECT,
         Diff::REFLECT, Diff::REFLECT};
    auto diff_solver = std::make_shared<Diff>(
        assembly, dz, bcs);

    int Nx = assembly->num_pins_x();
    int Ny = assembly->num_pins_y();
    auto N = Nx * Ny * dz.size();
    std::vector<double> temperatures(N, 1200.0);
    std::vector<double> densities(N, 0.7);
    std::vector<double> powers(N, 0.7);

    diff_solver->solve(temperatures, densities, powers);

    auto ref = powers[0];
    auto Nxy = Nx * Ny;
    for (int ij = 0; ij < Nxy; ++ij)
        EXPECT_SOFTEQ(ref, powers[ij], 1e-4);

    for (int ij = 0; ij < Nxy; ++ij)
        EXPECT_SOFTEQ(4*ref, powers[ij+Nxy], 1e-4);

    // Check normalization
    double nrm = 0.0;
    for (auto val : powers)
        nrm += val * val;
    nrm = std::sqrt(nrm);
    EXPECT_SOFTEQ(1.0, nrm, 1e-8);

    std::cout << "Power: ";
    for (auto val : powers)
        std::cout << val << " ";
    std::cout << std::endl;
}

TEST_F(TwoGroupDiffusionTest, singlepin_symmetry)
{
    std::vector<double> dz(20, 18.0);
    auto x_edges = {0.0, 1.26};
    auto y_edges = {0.0, 1.26};

    double height = 360.0;

    std::vector<AM::PIN_TYPE> pins =
        {AM::FUEL};
    auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);

    std::vector<Diff::BC_TYPE> bcs =
        {Diff::REFLECT, Diff::REFLECT,
         Diff::REFLECT, Diff::REFLECT,
         Diff::VACUUM,  Diff::VACUUM};
    auto diff_solver = std::make_shared<Diff>(
        assembly, dz, bcs);

    std::vector<double> temperatures(20, 1200.0);
    std::vector<double> densities(20, 0.7);
    std::vector<double> powers(20, 0.0);

    diff_solver->solve(temperatures, densities, powers);

    // Check symmetry
    for (size_t cell = 0; cell < dz.size(); ++cell)
        EXPECT_SOFTEQ(powers[cell], powers[dz.size()-1-cell], 1e-4);

    std::cout << "Power: ";
    for (auto val : powers)
        std::cout << val << " ";
    std::cout << std::endl;

    // Check normalization
    double nrm = 0.0;
    for (auto val : powers)
        nrm += val * val;
    nrm = std::sqrt(nrm);
    EXPECT_SOFTEQ(1.0, nrm, 1e-8);
}

TEST_F(TwoGroupDiffusionTest, singlepin_multitemp)
{
    auto x_edges = {0.0, 1.26};
    auto y_edges = {0.0, 1.26};
    std::vector<double> dz(30, 12.0);

    std::vector<AM::PIN_TYPE> pins =
        {AM::FUEL};
    double height = 360.0;
    auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);

    std::vector<Diff::BC_TYPE> bcs =
        {Diff::REFLECT, Diff::REFLECT,
         Diff::REFLECT, Diff::REFLECT,
         Diff::VACUUM,  Diff::VACUUM};
    auto diff_solver = std::make_shared<Diff>(
        assembly, dz, bcs);

    auto Nz = dz.size();
    std::vector<double> temperatures(Nz, 1000);
    std::vector<double> densities(Nz, 0.7);
    std::vector<double> powers(Nz, 0.0);

    for (int cell = 0; cell < Nz; ++cell)
    {
        double frac = static_cast<double>(cell) /
                      static_cast<double>(Nz);
        temperatures[cell] = 1000.0 + 500.0 * std::sin(3.14159 * frac);
        densities[cell] = 0.75 - frac * 0.1;
    }

    diff_solver->solve(temperatures, densities, powers);

    std::cout << "Power: ";
    for (auto val : powers)
        std::cout << val << " ";
    std::cout << std::endl;

    // Check normalization
    double nrm = 0.0;
    for (auto val : powers)
        nrm += val * val;
    nrm = std::sqrt(nrm);
    EXPECT_SOFTEQ(1.0, nrm, 1e-8);

    // Compare against heuristic
    std::vector<double> ref =
        {0.0506247, 0.145302,  0.224282,  0.282708,  0.319362,  0.335792,
         0.33529,   0.321938,  0.299865,  0.272758,  0.243606,  0.214631,
         0.187327,  0.162584,  0.140818,  0.122108,  0.106303,  0.0931171,
         0.0821846, 0.0731064, 0.0654724, 0.0588766, 0.0529253, 0.0472455,
         0.0414949, 0.0353784, 0.0286714, 0.0212484, 0.0131159, 0.00444117};
    EXPECT_VEC_SOFTEQ(ref, powers, 2e-4);
}

TEST_F(TwoGroupDiffusionTest, assembly)
{
    SKIP_TEST("Performance evaluation only");

    std::vector<double> x_edges(18, 0.0);
    std::vector<double> y_edges(18, 0.0);
    for (int i = 0; i < 17; ++i)
    {
        x_edges[i+1] = x_edges[i] + 1.26;
        y_edges[i+1] = y_edges[i] + 1.26;
    }
    std::vector<double> dz(60, 6.0);

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

    std::vector<Diff::BC_TYPE> bcs =
        {Diff::REFLECT, Diff::REFLECT,
         Diff::REFLECT, Diff::REFLECT,
         Diff::VACUUM,  Diff::VACUUM};
    auto diff_solver = std::make_shared<Diff>(
        assembly, dz, bcs);

    auto Nx = assembly->num_pins_x();
    auto Ny = assembly->num_pins_y();
    auto N = Nx * Ny * dz.size();
    std::vector<double> temperatures(N, 1000);
    std::vector<double> densities(N, 0.7);
    std::vector<double> powers(N, 0.0);

    diff_solver->solve(temperatures, densities, powers);
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstTwoGroupDiffusion.cc
//---------------------------------------------------------------------------//
