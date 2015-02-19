//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialFactory.cc
 * \author Steven Hamilton
 * \brief  Construct polynomial coefficients.
 */
//---------------------------------------------------------------------------//

#include <chrono>

#include "PolynomialFactory.hh"

#include "AleaTypedefs.hh"

#include "NeumannPolynomial.hh"
#include "ChebyshevPolynomial.hh"
#include "GmresPolynomial.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Construct polynomial class.
 *
 * Behavior is controlled by the following entries on the "Polynomial"
 * sublist:
 *  - polynomial_order(int) : order of polynomial to be constructed
 *  - polynomial_type(string) : ("neumann"), "chebyshev", "gmres"
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<Polynomial>
PolynomialFactory::buildPolynomial(Teuchos::RCP<const MATRIX> A,
                                   Teuchos::RCP<Teuchos::ParameterList> pl )
{
    TEUCHOS_ASSERT( pl != Teuchos::null );

    // A can be null for some polynomial types
    // Defer checking A to individual classes

    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");
    TEUCHOS_ASSERT( poly_pl != Teuchos::null );

    std::string poly_type = poly_pl->get("polynomial_type","neumann");
    TEUCHOS_ASSERT( poly_type == "neumann"   ||
                    poly_type == "chebyshev" ||
                    poly_type == "gmres" );

    Teuchos::RCP<Polynomial> poly;
    if( poly_type == "neumann" )
    {
        poly = Teuchos::rcp( new NeumannPolynomial(A,pl) );
    }
    else if( poly_type == "chebyshev" )
    {
        poly = Teuchos::rcp( new ChebyshevPolynomial(A,pl) );
    }
    else if( poly_type == "gmres" )
    {
        poly = Teuchos::rcp( new GmresPolynomial(A,pl) );
    }

    return poly;
}

} // namespace alea

