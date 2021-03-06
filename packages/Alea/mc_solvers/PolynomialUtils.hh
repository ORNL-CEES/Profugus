//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/PolynomialUtils.hh
 * \author Steven Hamilton
 * \brief  Pure static class for various polynomial-related utilities.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_PolynomialUtils_hh
#define Alea_mc_solvers_PolynomialUtils_hh

#include <cmath>

#include "Teuchos_Assert.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "AleaTypedefs.hh"

namespace alea
{
//---------------------------------------------------------------------------//
/*!
 * \class PolynomialUtils
 * \brief Helper functions for polynomial construction.
 */
//---------------------------------------------------------------------------//
class PolynomialUtils
{
  public:

    //! Binomial coefficient (choose function)
    static unsigned long choose(unsigned long n, unsigned long k);

    //! Compute coefficients for specified order Chebyshev polynomial
    static Teuchos::ArrayRCP<const SCALAR> getChebyshevCoefficients(LO n);
};

}

#endif // Alea_mc_solvers_PolynomialUtils_hh

