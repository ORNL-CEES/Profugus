//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NeumannPolynomial.hh
 * \author Steven Hamilton
 * \brief  NeumannPolynomial class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_NeumannPolynomial_hh
#define mc_solvers_NeumannPolynomial_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "PolynomialBasis.hh"
#include "Polynomial.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class NeumannPolynomial
 * \brief Neumann series polynomial.
 */
//---------------------------------------------------------------------------//
class NeumannPolynomial : public Polynomial
{

  public:

    // Constructor
    NeumannPolynomial(Teuchos::RCP<const MATRIX> A,
               Teuchos::RCP<Teuchos::ParameterList> pl);
};

}

#endif // mc_solvers_NeumannPolynomial_hh

