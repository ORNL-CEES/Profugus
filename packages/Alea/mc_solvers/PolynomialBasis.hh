//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialBasis.hh
 * \author Steven Hamilton
 * \brief  Class for evaluating and converting between polynomial bases.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_PolynomialBasis_hh
#define mc_solvers_PolynomialBasis_hh

#include <string>

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class PolynomialBasis
 * \brief Provide conversion between polynomial bases.
 *
 * It is common to define a polynomial in terms of the standard (or monomial)
 * basis, i.e.
 * \f[
        p(x) = \sum_{n=0}^N c_n x^n.
 * \f]
 * For numerical stability (particularly when dealing with polynomials of
 * matrices), however, it is often preferable to work with an alternate
 * basis, i.e.
 * \f[
        p(x) = \sum_{n=0}^N d_n (\alpha + \beta x)^n.
 * \f]
 * This class handles the conversion of the polynomial coefficients from
 * one basis to another.
 */
//---------------------------------------------------------------------------//
class PolynomialBasis
{

  public:

    // Constructor based on basis type ("neumann", "power", or "arbitrary")
    PolynomialBasis(std::string type);

    // Evaluate basis
    Teuchos::ArrayRCP<SCALAR>
    operator() (const Teuchos::ArrayRCP<const SCALAR> x, LO k) const;

    // Set basis coefficients
    void setBasisCoefficients(SCALAR alpha, SCALAR beta);

    //  Get basis coefficients \f$\alpha\f$ and \f$\beta\f$
    void getBasisCoefficients(SCALAR &alpha, SCALAR &beta) const;

    // Transform between bases
    Teuchos::ArrayRCP<SCALAR> transformBasis(
        const Teuchos::ArrayRCP<const SCALAR> other_coeffs,
        const PolynomialBasis &other_basis ) const;


  protected:

    // Get triangular matrix corresponding to transform to monomial basis
    Teuchos::SerialDenseMatrix<LO,SCALAR>
    getMonomialBasisTransform(LO k) const;

    enum TYPE {POWER, NEUMANN, ARBITRARY };

    TYPE   d_type;
    SCALAR d_alpha, d_beta;

};

}

#endif // mc_solvers_PolynomialBasis_hh

