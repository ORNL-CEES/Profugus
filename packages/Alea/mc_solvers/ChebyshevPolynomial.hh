//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/ChebyshevPolynomial.hh
 * \author Steven Hamilton
 * \brief  ChebyshevPolynomial class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_ChebyshevPolynomial_hh
#define Alea_mc_solvers_ChebyshevPolynomial_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "PolynomialBasis.hh"
#include "Polynomial.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class ChebyshevPolynomial
 * \brief Chebyshev series polynomial.
 */
//---------------------------------------------------------------------------//
class ChebyshevPolynomial : public Polynomial
{

  public:

    // Constructor
    ChebyshevPolynomial(Teuchos::RCP<const MATRIX> A,
               Teuchos::RCP<Teuchos::ParameterList> pl);

    //! \brief Access minimum eigenvalue.
    SCALAR getLambdaMin() const
    {
        REQUIRE( !SCALAR_TRAITS::isnaninf(d_lambda_min) );
        return d_lambda_min;
    }

    //! \brief Access maximum eigenvalue.
    SCALAR getLambdaMax() const
    {
        REQUIRE( !SCALAR_TRAITS::isnaninf(d_lambda_max) );
        return d_lambda_max;
    }

  private:

    void computeEigenvalues();

    // Eigenvalue data
    SCALAR d_lambda_min;
    SCALAR d_lambda_max;

};

}

#endif // Alea_mc_solvers_ChebyshevPolynomial_hh

