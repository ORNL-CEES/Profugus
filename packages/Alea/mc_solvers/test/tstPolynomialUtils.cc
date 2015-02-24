//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testPolynomialUtils
 * \author Steven Hamilton
 * \brief  Test of PolynomialUtils class.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "gtest/utils_gtest.hh"

#include "../AleaTypedefs.hh"
#include "../PolynomialUtils.hh"

namespace
{
unsigned long factorial(unsigned long n)
{
    unsigned long f = 1UL;
    for( unsigned long k=2UL; k<=n; ++k )
    {
        f *= k;
    }
    return f;
}
}

using namespace alea;

TEST(Choose, Basic)
{
    // Test choose function
    LO n, k;

    // Test with small values against naive definition
    for( n=0; n<10; ++n )
    {
        for( k=0; k<n; ++k )
        {
            unsigned long num = factorial(n);
            unsigned long denom1 = factorial(n-k);
            unsigned long denom2 = factorial(k);
            unsigned long expected = num / (denom1*denom2);
            EXPECT_EQ( expected, alea::PolynomialUtils::choose(n,k) );
        }
    }

    // Now test a few "large" values that will cause overflow
    //  with the naive implementation

    // choose(24,18) = 24!/(18! * 6!) = 24*23*22*21*20*19/6! = 23*22*7*2*19
    n = 24;
    k = 18;
    unsigned long expected = 23*22*7*2*19;
    EXPECT_EQ( expected, alea::PolynomialUtils::choose(n,k) );

    n = 32;
    k = 10;
    expected = 64512240;
    EXPECT_EQ( expected, alea::PolynomialUtils::choose(n,k) );
}

TEST(Chebyshev, Basic)
{
    // Test against a few specific cases
    LO n;
    Teuchos::ArrayRCP<const SCALAR> coeffs;


    // Expected values taken from wikipedia.org/wiki/Chebyshev_polynomials

    // n=1
    n=1;
    coeffs = alea::PolynomialUtils::getChebyshevCoefficients(n);
    EXPECT_EQ( n+1, coeffs.size() );
    EXPECT_DOUBLE_EQ( 0.0, coeffs[0] );
    EXPECT_DOUBLE_EQ( 1.0, coeffs[1] );

    // n=2
    n=2;
    coeffs = alea::PolynomialUtils::getChebyshevCoefficients(n);
    EXPECT_EQ( n+1, coeffs.size() );
    EXPECT_DOUBLE_EQ( -1.0, coeffs[0] );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[1] );
    EXPECT_DOUBLE_EQ(  2.0, coeffs[2] );

    // n=3
    n=3;
    coeffs = alea::PolynomialUtils::getChebyshevCoefficients(n);
    EXPECT_EQ( n+1, coeffs.size() );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[0] );
    EXPECT_DOUBLE_EQ( -3.0, coeffs[1] );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[2] );
    EXPECT_DOUBLE_EQ(  4.0, coeffs[3] );

    // n=4
    n=4;
    coeffs = alea::PolynomialUtils::getChebyshevCoefficients(n);
    EXPECT_EQ( n+1, coeffs.size() );
    EXPECT_DOUBLE_EQ(  1.0, coeffs[0] );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[1] );
    EXPECT_DOUBLE_EQ( -8.0, coeffs[2] );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[3] );
    EXPECT_DOUBLE_EQ(  8.0, coeffs[4] );

    // n=5
    n=5;
    coeffs = alea::PolynomialUtils::getChebyshevCoefficients(n);
    EXPECT_EQ( n+1, coeffs.size() );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[0] );
    EXPECT_DOUBLE_EQ(  5.0, coeffs[1] );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[2] );
    EXPECT_DOUBLE_EQ(-20.0, coeffs[3] );
    EXPECT_DOUBLE_EQ(  0.0, coeffs[4] );
    EXPECT_DOUBLE_EQ( 16.0, coeffs[5] );

    // n=8
    n=8;
    coeffs = alea::PolynomialUtils::getChebyshevCoefficients(n);
    EXPECT_EQ( n+1, coeffs.size() );
    EXPECT_DOUBLE_EQ(   1.0, coeffs[0] );
    EXPECT_DOUBLE_EQ(   0.0, coeffs[1] );
    EXPECT_DOUBLE_EQ( -32.0, coeffs[2] );
    EXPECT_DOUBLE_EQ(   0.0, coeffs[3] );
    EXPECT_DOUBLE_EQ( 160.0, coeffs[4] );
    EXPECT_DOUBLE_EQ(   0.0, coeffs[5] );
    EXPECT_DOUBLE_EQ(-256.0, coeffs[6] );
    EXPECT_DOUBLE_EQ(   0.0, coeffs[7] );
    EXPECT_DOUBLE_EQ( 128.0, coeffs[8] );
}

