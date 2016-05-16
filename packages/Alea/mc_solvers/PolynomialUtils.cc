//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/PolynomialUtils.cc
 * \author Steven Hamilton
 * \brief  Construct data necessary for Monte Carlo solvers.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iterator>

#include "PolynomialUtils.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Compute binomial coefficient
 *
 * \f$ \text{choose}(n,k) = \frac{n\!}{(n-k)\!k\!} \f$
 *
 * This algorithm is designed to reduce overflow issues relative to
 * a naive implementation from definition of the coefficient.
 */
//---------------------------------------------------------------------------//
unsigned long PolynomialUtils::choose(unsigned long n, unsigned long k)
{
    if (k > n) {
        return 0UL;
    }
    long r = 1UL;
    for (unsigned long d = 1UL; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute coefficients of Chebyshev polynomial of degree n.
 *
 * The Chebyshev polynomial of degree n can be written as
 * \f$ Tn(x) = sum_{k=0}^{floor(n/2)} d_k^{(n)} x^{n-2k} \f$
 * where
 * \f$ d_k^{(n)} = (-1)^k 2^{n-2k-1} \frac{n}{n-k}
 *   \left( \begin{array}{c} n-k \\ k \end{array} \right) \f$
 */
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<const SCALAR> PolynomialUtils::getChebyshevCoefficients(LO n)
{
    Teuchos::ArrayRCP<SCALAR> coeffs(n+1);

    for( LO k=0; k<=n/2; ++k )
    {
        // Compute coefficient
        SCALAR d = SCALAR_TRAITS::pow(-1.0,static_cast<SCALAR>(k));
        d *= SCALAR_TRAITS::pow(2.0,static_cast<SCALAR>(n-2*k-1));
        d *= static_cast<SCALAR>(n)/static_cast<SCALAR>(n-k);
        d *= static_cast<SCALAR>( PolynomialUtils::choose(n-k,k) );

        // This is coefficient of (n-2k) power of x
        coeffs[n-2*k] = d;
    }

    return coeffs;
}


} // namespace alea

