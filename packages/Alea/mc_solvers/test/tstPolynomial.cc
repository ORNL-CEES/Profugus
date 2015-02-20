//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testPolynomial
 * \author Steven Hamilton
 * \brief  Test of Polynomial class.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "gtest/utils_gtest.hh"

#include "../LinearSystem.hh"
#include "../LinearSystemFactory.hh"
#include "../Polynomial.hh"
#include "../PolynomialBasis.hh"
#include "../PolynomialFactory.hh"
#include "../NeumannPolynomial.hh"
#include "../ChebyshevPolynomial.hh"
#include "../GmresPolynomial.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace alea;

TEST(Polynomial, Basic)
{
    // Set problem parameters
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");
    poly_pl->set("reproducible_random",true);
    poly_pl->set("random_seed",54321);

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
    {
        std::vector<double> exp_vec = {3.0, -3.0, 1.0};
        for( int i=0; i<expected_size; ++i )
            EXPECT_DOUBLE_EQ(exp_vec[i],coeffs[i]);
    }

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
    EXPECT_EQ(expected_size,coeffs.size());
    {
        std::vector<double> exp_vec =
            {70.0/33.0, -48.0/33.0, 32.0/99.0};
        for( int i=0; i<expected_size; ++i )
            EXPECT_DOUBLE_EQ(exp_vec[i],coeffs[i]);
    }

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
    EXPECT_SOFTEQ( expected, lambda_min, 1e-6);
    SCALAR lambda_max = cheby_poly->getLambdaMax();
    expected = 1.99802681;
    EXPECT_SOFTEQ( expected, lambda_max, 1e-6);

    //
    // GMRES coefficients in power basis
    //

    std::cout << "GMRES coefficients in power basis" << std::endl;

    poly_pl->set("polynomial_type","gmres");
    poly_pl->set("polynomial_order",4);
    poly_pl->set("gmres_type","fom");

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
    {
        std::vector<double> exp_vec =
            { 120.3590920202, -588.4594349692, 915.9793178951,
             -562.3764435696,  118.7150058831};
        double tol = 1e-12;
        for( int i=0; i<expected_size; ++i )
            EXPECT_SOFTEQ(exp_vec[i],coeffs[i],tol);
    }

    std::cout << "GMRES power coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;

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

    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );
    // Compuare against MATLAB-computed values
    {
        std::vector<double> exp_vec =
            {317.3737917360, -6113.3991897250, 47858.3199791513,
             -194604.3988303378, 462323.7282052612, -684480.9307913361,
             650672.6102197397, -397405.5220882442, 150778.4680632685,
             -32345.3651561241, 2997.9784441415};
        double tol = 1e-12;
        for( int i=0; i<expected_size; ++i )
            EXPECT_SOFTEQ(exp_vec[i],coeffs[i],tol);
    }

    std::cout << "FOM coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;

    //
    // GMRES coefficients computed from normal equations solution
    //

    poly_pl->set("gmres_type","normal");

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
    {
        std::vector<double> exp_vec =
            {25.9658833933, -270.5756373551, 1466.4374282468,
             -4625.0162185862, 9110.8784873964, -11675.6152179777,
             9892.8402157345, -5499.6550235381, 1929.0035351178,
             -387.1279225560, 33.8812574180};
        double tol = 1e-2;
        for( int i=0; i<expected_size; ++i )
            EXPECT_SOFTEQ(exp_vec[i],coeffs[i],tol);
    }

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
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );

    {
        std::vector<double> exp_vec =
            {25.9715826978, -270.6862126003, 1467.2872852538,
             -4628.3926873274, 9118.7195336468, -11686.9826673866,
             9903.4421648398, -5506.0197798736, 1931.3811134933,
             -387.6308375889, 33.9272789271};
        double tol = 1e-9;
        for( int i=0; i<expected_size; ++i )
            EXPECT_SOFTEQ(exp_vec[i],coeffs[i],tol);
    }

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
    poly = PolynomialFactory::buildPolynomial(A,pl);
    EXPECT_TRUE( poly != Teuchos::null );

    // Make sure correct object was built
    gmres_poly = Teuchos::rcp_dynamic_cast<GmresPolynomial>(poly);
    TEUCHOS_ASSERT( gmres_poly != Teuchos::null );

    coeffs = poly->getCoeffs(power_basis);
    expected_size = 11;
    EXPECT_EQ( expected_size, coeffs.size() );

    {
        std::vector<double> exp_vec =
            {25.9715826978, -270.6862125900, 1467.2872851571,
             -4628.3926869256, 9118.7195327064, -11686.9826660298,
             9903.4421635866, -5506.0197791297, 1931.3811132186,
             -387.6308375315, 33.9272789218};

        double tol = 1e-10;
        for( int i=0; i<expected_size; ++i )
            EXPECT_SOFTEQ(exp_vec[i],coeffs[i],tol);
    }

    std::cout << "GMRES/Arnoldi coefficients: " << std::endl;
    for( int i=0; i<coeffs.size(); ++i )
        std::cout << i << " " << coeffs[i] << std::endl;
}

