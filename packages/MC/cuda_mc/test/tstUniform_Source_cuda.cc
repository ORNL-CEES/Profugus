//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstUniform_Source_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 15:22:15 2016
 * \brief  Test for Uniform_Source
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class Uniform_SourceTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS

  protected:
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Uniform_SourceTest, coincident)
{
    std::vector<double> geom_bounds = {-4.0, 2.0, -1.0, 1.0, 1.0, 12.0};

    double xlo =  0.0;
    double xhi =  1.0;
    double ylo = -1.0;
    double yhi =  1.0;
    double zlo =  3.0;
    double zhi =  5.0;
    std::vector<double> src_bounds = {xlo, xhi, ylo, yhi, zlo, zhi};

    int Np = 1024;
    std::vector<double> x_locs(Np);
    std::vector<double> y_locs(Np);
    std::vector<double> z_locs(Np);
    std::vector<double> mu_vals(Np);
    std::vector<double> eta_vals(Np);
    std::vector<double> xi_vals(Np);

    cuda_mc::Uniform_Source_Tester::test_source(
        geom_bounds,src_bounds,x_locs,y_locs,z_locs,
        mu_vals,eta_vals,xi_vals);

    // Now compute mean of position in each direction
    double x_expected = 0.5 * (xlo + xhi);
    double x_mean = 0.0;
    for( auto x : x_locs )
    {
        x_mean += x;
    }
    x_mean /= static_cast<double>(Np);

    double y_expected = 0.5 * (ylo + yhi);
    double y_mean = 0.0;
    for( auto y : y_locs )
    {
        y_mean += y;
    }
    y_mean /= static_cast<double>(Np);

    double z_expected = 0.5 * (zlo + zhi);
    double z_mean = 0.0;
    for( auto z : z_locs )
    {
        z_mean += z;
    }
    z_mean /= static_cast<double>(Np);

    std::cout << "Expecting mean position of (" << x_expected
        << ", " << y_expected << ", " << z_expected << ")" << std::endl;
    std::cout << "Actual mean x = " << x_mean << std::endl;
    std::cout << "Actual mean y = " << y_mean << std::endl;
    std::cout << "Actual mean z = " << z_mean << std::endl;

    double sqrt_N = std::sqrt(static_cast<double>(Np));
    EXPECT_TRUE( std::abs(x_mean - x_expected) < (xhi - xlo)/sqrt_N );
    EXPECT_TRUE( std::abs(y_mean - y_expected) < (yhi - ylo)/sqrt_N );
    EXPECT_TRUE( std::abs(z_mean - z_expected) < (zhi - zlo)/sqrt_N );

    // Now compute mean of each direction component
    double mu_mean = 0.0;
    for( auto mu : mu_vals)
    {
        mu_mean += mu;
    }
    mu_mean /= static_cast<double>(Np);

    double eta_mean = 0.0;
    for( auto eta : eta_vals )
    {
        eta_mean += eta;
    }
    eta_mean /= static_cast<double>(Np);

    double xi_mean = 0.0;
    for( auto xi : xi_vals )
    {
        xi_mean += xi;
    }
    xi_mean /= static_cast<double>(Np);

    std::cout << "Expecting mean angle of (0,0,0)" << std::endl;
    std::cout << "Actual mean mu = " << mu_mean << std::endl;
    std::cout << "Actual mean eta = " << eta_mean << std::endl;
    std::cout << "Actual mean xi = " << xi_mean << std::endl;

    // All angle sampling is isotropic so expected mean is zero
    double mean_expected = 0.0;
    EXPECT_TRUE( std::abs(mu_mean  - mean_expected) < 1.0/sqrt_N );
    EXPECT_TRUE( std::abs(eta_mean - mean_expected) < 1.0/sqrt_N );
    EXPECT_TRUE( std::abs(xi_mean  - mean_expected) < 1.0/sqrt_N );
}

TEST_F(Uniform_SourceTest, host)
{
    std::vector<double> geom_bounds = {-4.0, 2.0, -1.0, 1.0, 1.0, 12.0};

    double xlo =  0.0;
    double xhi =  1.0;
    double ylo = -1.0;
    double yhi =  1.0;
    double zlo =  3.0;
    double zhi =  5.0;
    std::vector<double> src_bounds = {xlo, xhi, ylo, yhi, zlo, zhi};

    int Np = 1024;
    std::vector<double> x_locs(Np);
    std::vector<double> y_locs(Np);
    std::vector<double> z_locs(Np);
    std::vector<double> mu_vals(Np);
    std::vector<double> eta_vals(Np);
    std::vector<double> xi_vals(Np);

    cuda_mc::Uniform_Source_Tester::test_host_api(
        geom_bounds,src_bounds,x_locs,y_locs,z_locs,
        mu_vals,eta_vals,xi_vals);


    // Now compute mean of position in each direction
    double x_expected = 0.5 * (xlo + xhi);
    double x_mean = 0.0;
    for( auto x : x_locs )
    {
        x_mean += x;
    }
    x_mean /= static_cast<double>(Np);

    double y_expected = 0.5 * (ylo + yhi);
    double y_mean = 0.0;
    for( auto y : y_locs )
    {
        y_mean += y;
    }
    y_mean /= static_cast<double>(Np);

    double z_expected = 0.5 * (zlo + zhi);
    double z_mean = 0.0;
    for( auto z : z_locs )
    {
        z_mean += z;
    }
    z_mean /= static_cast<double>(Np);

    std::cout << "Expecting mean position of (" << x_expected
        << ", " << y_expected << ", " << z_expected << ")" << std::endl;
    std::cout << "Actual mean x = " << x_mean << std::endl;
    std::cout << "Actual mean y = " << y_mean << std::endl;
    std::cout << "Actual mean z = " << z_mean << std::endl;

    double sqrt_N = std::sqrt(static_cast<double>(Np));
    EXPECT_TRUE( std::abs(x_mean - x_expected) < (xhi - xlo)/sqrt_N );
    EXPECT_TRUE( std::abs(y_mean - y_expected) < (yhi - ylo)/sqrt_N );
    EXPECT_TRUE( std::abs(z_mean - z_expected) < (zhi - zlo)/sqrt_N );

    // Now compute mean of each direction component
    double mu_mean = 0.0;
    for( auto mu : mu_vals)
    {
        mu_mean += mu;
    }
    mu_mean /= static_cast<double>(Np);

    double eta_mean = 0.0;
    for( auto eta : eta_vals )
    {
        eta_mean += eta;
    }
    eta_mean /= static_cast<double>(Np);

    double xi_mean = 0.0;
    for( auto xi : xi_vals )
    {
        xi_mean += xi;
    }
    xi_mean /= static_cast<double>(Np);

    std::cout << "Expecting mean angle of (0,0,0)" << std::endl;
    std::cout << "Actual mean mu = " << mu_mean << std::endl;
    std::cout << "Actual mean eta = " << eta_mean << std::endl;
    std::cout << "Actual mean xi = " << xi_mean << std::endl;

    // All angle sampling is isotropic so expected mean is zero
    double mean_expected = 0.0;
    EXPECT_TRUE( std::abs(mu_mean  - mean_expected) < 1.0/sqrt_N );
    EXPECT_TRUE( std::abs(eta_mean - mean_expected) < 1.0/sqrt_N );
    EXPECT_TRUE( std::abs(xi_mean  - mean_expected) < 1.0/sqrt_N );
}

//---------------------------------------------------------------------------//
//                 end of tstUniform_Source.cc
//---------------------------------------------------------------------------//
