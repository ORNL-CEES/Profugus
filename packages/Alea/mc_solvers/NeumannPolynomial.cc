//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/NeumannPolynomial.cc
 * \author Steven Hamilton
 * \brief  NeumannPolynomial class definitions.
 */
//---------------------------------------------------------------------------//

#include "NeumannPolynomial.hh"
#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Construct Neumann polynomial.
 */
//---------------------------------------------------------------------------//
NeumannPolynomial::NeumannPolynomial(Teuchos::RCP<const MATRIX> A,
                                     Teuchos::RCP<Teuchos::ParameterList> pl)
  : Polynomial(A,pl)
{
    // Damped Neumann polynomial basis is given by (I - omega*A)
    SCALAR damp = b_poly_pl->get("neumann_damping",1.0);

    // Create native basis depending on damping parameter
    if( damp == 1.0 )
    {
        if( b_verbosity >= LOW )
            std::cout << "Creating Neumann polynomial coefficients"
                << " of order " << b_m << std::endl;

        b_native_basis = Teuchos::rcp( new PolynomialBasis("neumann") );
    }
    else
    {
        if( b_verbosity >= LOW )
            std::cout << "Creating damped Neumann polynomial coefficients"
                << " of order " << b_m << std::endl;

        b_native_basis = Teuchos::rcp( new PolynomialBasis("arbitrary") );
        b_native_basis->setBasisCoefficients(1.0,-damp);
    }
    CHECK( b_native_basis != Teuchos::null );

    b_coeffs.resize(b_m+1);
    std::fill( b_coeffs.begin(), b_coeffs.end(), damp );
}

} // namespace alea

