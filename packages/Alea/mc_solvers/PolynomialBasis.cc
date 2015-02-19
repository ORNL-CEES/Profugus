//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialBasis.cc
 * \author Steven Hamilton
 * \brief  Construct data necessary for Monte Carlo solvers.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <limits>
#include <vector>

#include "PolynomialBasis.hh"

#include "Teuchos_Assert.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param type Type of basis corresponding to this object.
 *
 * Accepted values are:
 *  - "power" : The monomial basis \f$x^n\f$
 *  - "neumann" : The Neumann basis \f$(1-x)^n\f$
 *  - "arbitrary" : An arbitrary basis \f$(\alpha + \beta x)^n\f$
 *
 * If the basis type is set to "arbitrary", the values of \f$\alpha\f$ and
 * \f$\beta\f$ MUST be set via a call to setBasisCoefficients before using
 * this object.
 */
//---------------------------------------------------------------------------//
PolynomialBasis::PolynomialBasis(std::string type)
{
    if( type == "power" )
    {
        d_alpha = 0.0;
        d_beta  = 1.0;
        d_type  = POWER;
    }
    else if( type == "neumann" )
    {
        d_alpha =  1.0;
        d_beta  = -1.0;
        d_type  = NEUMANN;
    }
    else if( type == "arbitrary" )
    {
        d_alpha = SCALAR_TRAITS::nan();
        d_beta  = SCALAR_TRAITS::nan();
        d_type  = ARBITRARY;
    }
    else
    {
        TEUCHOS_ASSERT(false);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute \f$(\alpha + \beta x)^k\f$ for each value of x
 *
 * \param x Array of values of x to be evaluated.
 * \param k Exponent of evaluation.
 */
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<SCALAR> PolynomialBasis::operator()(
        const Teuchos::ArrayRCP<const SCALAR> x, LO k) const
{
    Teuchos::ArrayRCP<SCALAR> y(x.size());
    for( LO i=0; i<x.size(); ++i )
    {
        SCALAR base = d_alpha + d_beta*x[i];
        y[i] = SCALAR_TRAITS::pow(base,static_cast<SCALAR>(k));
    }
    return y;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set basis coefficients \f$\alpha\f$ and \f$\beta\f$
 *
 * Only valid if polynomial basis type is "arbitrary"
 */
//---------------------------------------------------------------------------//
void PolynomialBasis::setBasisCoefficients(SCALAR alpha, SCALAR beta)
{
    TEUCHOS_ASSERT(d_type==ARBITRARY);
    d_alpha = alpha;
    d_beta  = beta;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get basis coefficients \f$\alpha\f$ and \f$\beta\f$
 */
//---------------------------------------------------------------------------//
void PolynomialBasis::getBasisCoefficients(SCALAR &alpha, SCALAR &beta) const
{
    alpha = d_alpha;
    beta  = d_beta;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transform coefficients from another basis to this one.
 *
 * \param other_coeffs Array of coefficients with respect to other_basis
 * \param other_basis  PolynomialBasis object corresponding to other basis.
 *
 * Given a set of coefficients and a basis for those coefficients, compute
 * a new set of coefficients in the current basis.  This transformation is
 * computed by first converting the coefficients from the other basis to
 * the standard (monomial) basis via an upper triangular matrix-vector
 * product and then converting from the standard basis to the new basis
 * via an upper triangular linear solve, i.e.
 * \f$ c_n = M_n^{-1} M_o c_o \f$,
 * where \f$M_o\f$ and \f$M_n\f$ are the matrices corresponding to tranformations
 * from the other and new bases, respectively, to the standard basis,
 * and \f$c_o\f$ and \f$c_n\f$ are the polynomial coefficients with respect
 * to the corresponding bases.
 */
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<SCALAR> PolynomialBasis::transformBasis(
        const Teuchos::ArrayRCP<const SCALAR> other_coeffs,
        const PolynomialBasis &other_basis ) const
{
    TEUCHOS_ASSERT( d_alpha != SCALAR_TRAITS::nan() );
    TEUCHOS_ASSERT( d_beta  != SCALAR_TRAITS::nan() );

    SCALAR other_alpha, other_beta;
    other_basis.getBasisCoefficients(other_alpha,other_beta);
    TEUCHOS_ASSERT( other_alpha != SCALAR_TRAITS::nan() );
    TEUCHOS_ASSERT( other_beta  != SCALAR_TRAITS::nan() );

    // If this basis is equal to other basis, coefficients are unchanged
    if( (other_alpha==d_alpha) && (other_beta==d_beta) )
    {
        // For const-correctness, need to do a copy
        Teuchos::ArrayRCP<SCALAR> copy;
        copy.deepCopy( other_coeffs() );
        return copy;
    }

    LO N = other_coeffs.size();

    // Get triangular basis transformation matrices
    Teuchos::SerialDenseMatrix<LO,SCALAR> Ms =
        other_basis.getMonomialBasisTransform(N);
    Teuchos::SerialDenseMatrix<LO,SCALAR> Mt =
        this->getMonomialBasisTransform(N);

    Teuchos::BLAS<LO,SCALAR> blas;
    Teuchos::ArrayRCP<SCALAR> new_coeffs(N);
    new_coeffs.deepCopy(other_coeffs());

    // First tmp = M1 * other_coeffs
    blas.TRMV(Teuchos::UPPER_TRI,Teuchos::NO_TRANS,Teuchos::NON_UNIT_DIAG,
              N,Ms.values(),Ms.stride(),new_coeffs.get(),1);

    // Next new_coeffs = M2 \ tmp
    blas.TRSM(Teuchos::LEFT_SIDE,Teuchos::UPPER_TRI,Teuchos::NO_TRANS,
              Teuchos::NON_UNIT_DIAG,N,1,1.0,Mt.values(),Mt.stride(),
              new_coeffs.get(),N);

    return new_coeffs;
}

//---------------------------------------------------------------------------//
// PROTECTED MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Get matrix corresponding to transformation to standard basis
//---------------------------------------------------------------------------//
Teuchos::SerialDenseMatrix<LO,SCALAR>
PolynomialBasis::getMonomialBasisTransform(LO k) const
{
    bool zero_out = true;
    Teuchos::SerialDenseMatrix<LO,SCALAR> T(k,k,zero_out);
    T(0,0) = SCALAR_TRAITS::one();
    for( LO i=1; i<k; ++i )
    {
        for( LO j=0; j<i; ++j )
        {
            T(j,i)   += d_alpha * T(j,i-1);
            T(j+1,i) += d_beta * T(j,i-1);
        }
    }

    return T;
}

} // namespace alea

