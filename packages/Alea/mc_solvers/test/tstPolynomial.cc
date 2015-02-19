//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testPolynomial
 * \author Steven Hamilton
 * \brief  Test of Polynomial class.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "gtest/utils_gtest.hh"

#include "LinearSystem.hh"
#include "LinearSystemFactory.hh"
#include "Polynomial.hh"
#include "PolynomialBasis.hh"
#include "PolynomialFactory.hh"
#include "NeumannPolynomial.hh"
#include "ChebyshevPolynomial.hh"
#include "GmresPolynomial.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

using namespace alea;

TEST(Polynomial, Basic)

    // Set problem parameters
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");

    // Construct Map
    LO N = 50;
    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",N);
    mat_pl->set("scaling_type","diagonal");

    // Build linear system and get matrix
    Teuchos::RCP<alea::LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();

    //
    // Test standard Neumann coefficients
    //

    std::cout << "Neumann coefficient in Neumann basis" << std::endl;

    poly_pl->set("polynomial_type","neumann");
    poly_pl->set("polynomial_order",10);

    Teuchos::RCP<Polynomial> poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    Teuchos::RCP<NeumannPolynomial> neumann_poly =
        Teuchos::rcp_dynamic_cast<NeumannPolynomial>(poly);
    TEUCHOS_ASSERT( neumann_poly != Teuchos::null );

    PolynomialBasis neumann_basis("neumann");
    Teuchos::ArrayRCP<const SCALAR> coeffs = poly->getCoeffs(neumann_basis);
    LO expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );
    SCALAR expected = 1.0;
    for( size_t i=0; i<11; ++i )
    {
        EXPECT_DOUBLE_EQ( expected, coeffs[i] );
    }

    //
    // Neumann coefficients in power basis
    //

    std::cout << "Neumann coefficient in power basis" << std::endl;

    poly_pl->set("polynomial_type","neumann");
    poly_pl->set("polynomial_order",2);

    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    neumann_poly = Teuchos::rcp_dynamic_cast<NeumannPolynomial>(poly);
    TEUCHOS_ASSERT( neumann_poly != Teuchos::null );

    // Test coefficients against analytic values
    PolynomialBasis power_basis("power");
    coeffs = poly->getCoeffs(power_basis);
    expected_size = 3;
    EXPECT_EQ( expected_size, coeffs.size() );
    expected = 3.0;
    EXPECT_DOUBLE_EQ( expected, coeffs[0] );
    expected = -3.0;
    EXPECT_DOUBLE_EQ( expected, coeffs[1] );
    expected = 1.0;
    EXPECT_DOUBLE_EQ( expected, coeffs[2] );

    //
    // Chebyshev coefficients in power basis
    //

    std::cout << "Chebyshev coefficient in power basis" << std::endl;

    poly_pl->set("polynomial_type","chebyshev");
    poly_pl->set("polynomial_order",2);

    // Don't compute eigenvalues, just use prescribed values
    // These are not the true eigenvalues but we're just making
    // sure that we get the right coefficients for given eigenvalues
    poly_pl->set("compute_eigenvalues",false);
    poly_pl->set("lambda_min",1.0);
    poly_pl->set("lambda_max",2.0);

    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    Teuchos::RCP<ChebyshevPolynomial> cheby_poly =
        Teuchos::rcp_dynamic_cast<ChebyshevPolynomial>(poly);
    TEUCHOS_ASSERT( cheby_poly != Teuchos::null );

    // Test coefficients against analytic values
    coeffs = poly->getCoeffs(power_basis);
    expected_size = 3;
    EXPECT_EQ( expected_size, coeffs.size() );
    expected = 70.0/33.0;
    EXPECT_DOUBLE_EQ( expected, coeffs[0] );
    expected = -48.0/33.0;
    EXPECT_DOUBLE_EQ( expected, coeffs[1] );
    expected = 32.0/99.0;
    EXPECT_DOUBLE_EQ( expected, coeffs[2] );

    //
    // Check eigenvalue calculation in Chebyshev construction
    //

    std::cout << "Eigenvalue calculation" << std::endl;

    poly_pl->set("polynomial_type","chebyshev");
    poly_pl->set("polynomial_order",2);
    poly_pl->set("compute_eigenvalues",true);

    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    cheby_poly = Teuchos::rcp_dynamic_cast<ChebyshevPolynomial>(poly);
    TEUCHOS_ASSERT( cheby_poly != Teuchos::null );

    // Compuare against MATLAB-computed values
    SCALAR lambda_min = cheby_poly->getLambdaMin();
    expected = 1.9731936e-3;
    EXPECT_NEAR( expected, lambda_min, 1e-6);
    SCALAR lambda_max = cheby_poly->getLambdaMax();
    expected = 1.99802681;
    EXPECT_NEAR( expected, lambda_max, 1e-6);

    //
    // GMRES coefficients in power basis
    //

    std::cout << "GMRES coefficients in power basis" << std::endl;

    poly_pl->set("polynomial_type","gmres");
    poly_pl->set("polynomial_order",4);
    poly_pl->set("gmres_type","fom");

    // GMRES polynomial uses Teuchos::ScalarTraits to select
    // random values for initial Krylov vector
    // Set particular seed here for reproducibility in test
    SCALAR_TRAITS::seedrandom(10101);
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    Teuchos::RCP<GmresPolynomial> gmres_poly =
        Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 5;
    EXPECT_EQ( expected_size, coeffs.size() );
    // Compare against MATLAB-computed values
    expected = 36.2463717933919;
    EXPECT_NEAR( expected, coeffs[0], 1e-10 );
    expected = -156.7745506438811;
    EXPECT_NEAR( expected, coeffs[1], 1e-10 );
    expected = 225.4189277423314;
    EXPECT_NEAR( expected, coeffs[2], 1e-10 );
    expected = -130.0547044559239;
    EXPECT_NEAR( expected, coeffs[3], 1e-10 );
    expected = 26.0855267874441;
    EXPECT_NEAR( expected, coeffs[4], 1e-10 );

    //
    // GMRES coefficients for matrix with complex eigenvalues
    //

    std::cout << "GMRES coefficients with complex spectrum" << std::endl;

    mat_pl->set("matrix_type","convection-diffusion");
    mat_pl->set("matrix_size",100);
    mat_pl->set("viscosity",0.02);

    // Build linear system and get matrix
    system = alea::LinearSystemFactory::buildLinearSystem(pl);
    A = system->getMatrix();

    poly_pl->set("polynomial_type","gmres");
    poly_pl->set("polynomial_order",10);
    poly_pl->set("gmres_type","fom");

    // GMRES polynomial uses Teuchos::ScalarTraits to select
    // random values for initial Krylov vector
    // Set particular seed here for reproducibility in test
    SCALAR_TRAITS::seedrandom(10101);
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );
    // Compuare against MATLAB-computed values
    expected = 41.85651066835;
    EXPECT_NEAR( expected, coeffs[0], 1e-8 );
    expected = -615.66330986050;
    EXPECT_NEAR( expected, coeffs[1], 1e-8 );
    expected = 4136.94716550727;
    EXPECT_NEAR( expected, coeffs[2], 1e-8 );
    expected = -15242.13854812877;
    EXPECT_NEAR( expected, coeffs[3], 1e-8 );
    expected = 33939.59056325329;
    EXPECT_NEAR( expected, coeffs[4], 1e-8 );
    expected = -48059.76677361113;
    EXPECT_NEAR( expected, coeffs[5], 1e-8 );
    expected = 44223.10325663088;
    EXPECT_NEAR( expected, coeffs[6], 1e-8 );
    expected = -26333.51323537354;
    EXPECT_NEAR( expected, coeffs[7], 1e-8 );
    expected = 9784.04803444970;
    EXPECT_NEAR( expected, coeffs[8], 1e-8 );
    expected = -2061.15848172241;
    EXPECT_NEAR( expected, coeffs[9], 1e-8 );
    expected = 187.95599017072;
    EXPECT_NEAR( expected, coeffs[10], 1e-8 );

    std::cout << "FOM coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;

    //
    // GMRES coefficients computed from normal equations solution
    //

    poly_pl->set("gmres_type","normal");

    // GMRES polynomial uses Teuchos::ScalarTraits to select
    // random values for initial Krylov vector
    // Set particular seed here for reproducibility in test
    SCALAR_TRAITS::seedrandom(10101);
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );

    // Due to poor conditioning of linear system,
    //  only first few digits can be expected to match
    expected =  18.01626600611;
    EXPECT_NEAR( expected, coeffs[0], 1e-2*std::abs(expected) );
    expected = -159.07098510860;
    EXPECT_NEAR( expected, coeffs[1], 1e-2*std::abs(expected) );
    expected =  777.60948752050;
    EXPECT_NEAR( expected, coeffs[2], 1e-2*std::abs(expected) );
    expected = -2294.89032307120;
    EXPECT_NEAR( expected, coeffs[3], 1e-2*std::abs(expected) );
    expected =  4324.66928714378;
    EXPECT_NEAR( expected, coeffs[4], 1e-2*std::abs(expected) );
    expected = -5370.64743029258;
    EXPECT_NEAR( expected, coeffs[5], 1e-2*std::abs(expected) );
    expected =  4443.23962490253;
    EXPECT_NEAR( expected, coeffs[6], 1e-2*std::abs(expected) );
    expected = -2422.62831424262;
    EXPECT_NEAR( expected, coeffs[7], 1e-2*std::abs(expected) );
    expected =  835.66183653430;
    EXPECT_NEAR( expected, coeffs[8], 1e-2*std::abs(expected) );
    expected = -165.20776507683;
    EXPECT_NEAR( expected, coeffs[9], 1e-2*std::abs(expected) );
    expected =  14.25915717817;
    EXPECT_NEAR( expected, coeffs[10], 1e-2*std::abs(expected) );

    std::cout << "GMRES/NE coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;

    //
    // GMRES coefficients computed from QR factorization
    //

    poly_pl->set("gmres_type","qr");

    // GMRES polynomial uses Teuchos::ScalarTraits to select
    // random values for initial Krylov vector
    // Set particular seed here for reproducibility in test
    SCALAR_TRAITS::seedrandom(10101);
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );

    // Due to poor conditioning of linear system,
    //  only first few digits can be expected to match
    expected =  18.01626600611;
    EXPECT_NEAR( expected, coeffs[0], 1e-3*std::abs(expected) );
    expected = -159.07098510860;
    EXPECT_NEAR( expected, coeffs[1], 1e-3*std::abs(expected) );
    expected =  777.60948752050;
    EXPECT_NEAR( expected, coeffs[2], 1e-3*std::abs(expected) );
    expected = -2294.89032307120;
    EXPECT_NEAR( expected, coeffs[3], 1e-3*std::abs(expected) );
    expected =  4324.66928714378;
    EXPECT_NEAR( expected, coeffs[4], 1e-3*std::abs(expected) );
    expected = -5370.64743029258;
    EXPECT_NEAR( expected, coeffs[5], 1e-3*std::abs(expected) );
    expected =  4443.23962490253;
    EXPECT_NEAR( expected, coeffs[6], 1e-3*std::abs(expected) );
    expected = -2422.62831424262;
    EXPECT_NEAR( expected, coeffs[7], 1e-3*std::abs(expected) );
    expected =  835.66183653430;
    EXPECT_NEAR( expected, coeffs[8], 1e-3*std::abs(expected) );
    expected = -165.20776507683;
    EXPECT_NEAR( expected, coeffs[9], 1e-3*std::abs(expected) );
    expected =  14.25915717817;
    EXPECT_NEAR( expected, coeffs[10], 1e-3*std::abs(expected) );

    std::cout << "GMRES/QR coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;

    //
    // GMRES coefficients computed from Arnoldi
    //

    poly_pl->set("gmres_type","arnoldi");

    // GMRES polynomial uses Teuchos::ScalarTraits to select
    // random values for initial Krylov vector
    // Set particular seed here for reproducibility in test
    SCALAR_TRAITS::seedrandom(10101);
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );

    // Due to poor conditioning of linear system,
    //  only first few digits can be expected to match
    expected =  18.01626600611;
    EXPECT_NEAR( expected, coeffs[0], 1e-3*std::abs(expected) );
    expected = -159.07098510860;
    EXPECT_NEAR( expected, coeffs[1], 1e-3*std::abs(expected) );
    expected =  777.60948752050;
    EXPECT_NEAR( expected, coeffs[2], 1e-3*std::abs(expected) );
    expected = -2294.89032307120;
    EXPECT_NEAR( expected, coeffs[3], 1e-3*std::abs(expected) );
    expected =  4324.66928714378;
    EXPECT_NEAR( expected, coeffs[4], 1e-3*std::abs(expected) );
    expected = -5370.64743029258;
    EXPECT_NEAR( expected, coeffs[5], 1e-3*std::abs(expected) );
    expected =  4443.23962490253;
    EXPECT_NEAR( expected, coeffs[6], 1e-3*std::abs(expected) );
    expected = -2422.62831424262;
    EXPECT_NEAR( expected, coeffs[7], 1e-3*std::abs(expected) );
    expected =  835.66183653430;
    EXPECT_NEAR( expected, coeffs[8], 1e-3*std::abs(expected) );
    expected = -165.20776507683;
    EXPECT_NEAR( expected, coeffs[9], 1e-3*std::abs(expected) );
    expected =  14.25915717817;
    EXPECT_NEAR( expected, coeffs[10], 1e-3*std::abs(expected) );

    std::cout << "GMRES/Arnoldi coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;
}

