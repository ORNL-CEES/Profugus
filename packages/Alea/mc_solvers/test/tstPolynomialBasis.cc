//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testPolynomialBasis
 * \author Steven Hamilton
 * \brief  Test of PolynomialBasis class.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "gtest/utils_gtest.hh"

#include "../AleaTypedefs.hh"
#include "../PolynomialBasis.hh"

using namespace alea;

TEST(PolynomialBasis, Basic)
{
    LO num_values=5;
    Teuchos::ArrayRCP<SCALAR> x(num_values);
    for( LO i=0; i<num_values; ++i )
    {
        x[i] = static_cast<SCALAR>(i);
    }

    // Build and test power basis
    std::string type = "power";
    alea::PolynomialBasis power_basis(type);
    SCALAR alpha=0.0, beta=0.0;
    power_basis.getBasisCoefficients(alpha,beta);
    EXPECT_DOUBLE_EQ( 0.0, alpha );
    EXPECT_DOUBLE_EQ( 1.0, beta );

    Teuchos::ArrayRCP<const SCALAR> y;
    LO num_powers = 10;
    for( LO k=0; k<num_powers; ++k )
    {
        y = power_basis(x,k);
        for( LO i=0; i<num_values; ++i )
        {
            EXPECT_DOUBLE_EQ( SCALAR_TRAITS::pow(static_cast<SCALAR>(i),k), y[i] );
        }
    }

    // Build and test Neumann basis
    type = "neumann";
    alea::PolynomialBasis neumann_basis(type);
    alpha = 0.0;
    beta = 0.0;
    neumann_basis.getBasisCoefficients(alpha,beta);
    EXPECT_DOUBLE_EQ(  1.0, alpha );
    EXPECT_DOUBLE_EQ( -1.0, beta );

    for( LO k=0; k<num_powers; ++k )
    {
        y = neumann_basis(x,k);
        for( LO i=0; i<num_values; ++i )
        {
            EXPECT_DOUBLE_EQ( SCALAR_TRAITS::pow(1.0-static_cast<SCALAR>(i),k), y[i] );
        }
    }

    //
    // Test transformations
    //

    // Power -> Power
    Teuchos::ArrayRCP<SCALAR> input_coeffs(3);
    Teuchos::ArrayRCP<SCALAR> output_coeffs(3);
    input_coeffs[0] =  1.0;
    input_coeffs[1] = -2.0;
    input_coeffs[2] =  3.0;
    output_coeffs = power_basis.transformBasis(input_coeffs,power_basis);
    EXPECT_DOUBLE_EQ( input_coeffs[0], output_coeffs[0] );
    EXPECT_DOUBLE_EQ( input_coeffs[1], output_coeffs[1] );
    EXPECT_DOUBLE_EQ( input_coeffs[2], output_coeffs[2] );

    // Neumann -> Power
    output_coeffs = power_basis.transformBasis(input_coeffs,neumann_basis);
    EXPECT_DOUBLE_EQ(  2.0, output_coeffs[0] );
    EXPECT_DOUBLE_EQ( -4.0, output_coeffs[1] );
    EXPECT_DOUBLE_EQ(  3.0, output_coeffs[2] );

    // Power -> Neumann
    input_coeffs[0] =  2.0;
    input_coeffs[1] = -4.0;
    input_coeffs[2] =  3.0;
    output_coeffs = neumann_basis.transformBasis(input_coeffs,power_basis);
    EXPECT_DOUBLE_EQ(  1.0, output_coeffs[0] );
    EXPECT_DOUBLE_EQ( -2.0, output_coeffs[1] );
    EXPECT_DOUBLE_EQ(  3.0, output_coeffs[2] );

    // Arbitrary to Arbitrary
    type = "arbitrary";
    alea::PolynomialBasis arb1(type);
    arb1.setBasisCoefficients(1.0,3.0);
    alea::PolynomialBasis arb2(type);
    arb2.setBasisCoefficients(2.0,-1.0);

    input_coeffs.resize(4);
    output_coeffs.resize(4);
    input_coeffs[0] =  2.0;
    input_coeffs[1] = -3.0;
    input_coeffs[2] = -4.0;
    input_coeffs[3] =  5.0;
    output_coeffs = arb2.transformBasis(input_coeffs,arb1);
    EXPECT_DOUBLE_EQ(  1500.0, output_coeffs[0] );
    EXPECT_DOUBLE_EQ( -2028.0, output_coeffs[1] );
    EXPECT_DOUBLE_EQ(   909.0, output_coeffs[2] );
    EXPECT_DOUBLE_EQ(  -135.0, output_coeffs[3] );

    // And back
    input_coeffs[0] =  1500.0;
    input_coeffs[1] = -2028.0;
    input_coeffs[2] =   909.0;
    input_coeffs[3] =  -135.0;
    output_coeffs = arb1.transformBasis(input_coeffs,arb2);
    EXPECT_DOUBLE_EQ(  2.0, output_coeffs[0] );
    EXPECT_DOUBLE_EQ( -3.0, output_coeffs[1] );
    EXPECT_DOUBLE_EQ( -4.0, output_coeffs[2] );
    EXPECT_DOUBLE_EQ(  5.0, output_coeffs[3] );
}

