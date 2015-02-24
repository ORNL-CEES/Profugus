//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialPreconditioner.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "PolynomialPreconditioner.hh"
#include "PolynomialFactory.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * See documentation for Polynomial to determine polynomial options.
 */
//---------------------------------------------------------------------------//
PolynomialPreconditioner::PolynomialPreconditioner(
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl )
    : AleaSolver(A,pl)
{
    // Build polynomial and extract coefficients
    Teuchos::RCP<Polynomial> poly = PolynomialFactory::buildPolynomial(A,pl);
    REQUIRE( poly != Teuchos::null );
    d_coeffs = poly->getCoeffs();
    REQUIRE( d_coeffs != Teuchos::null );
    REQUIRE( d_coeffs.size() > 0 );
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Perform polynomial preconditioning.
 */
//---------------------------------------------------------------------------//
void PolynomialPreconditioner::applyImpl(const MV &x, MV &y) const
{
    // For now we only support operating on a single vector
    REQUIRE( x.getNumVectors() == 1 );

    bool init_to_zero = true;
    MV tmp_result(y.getMap(),y.getNumVectors(),init_to_zero);

    int num_coeffs = d_coeffs.size();
    y.update(d_coeffs[num_coeffs-1],x,0.0);

    // Use Horner's method to apply polynomial
    for( int i=num_coeffs-1; i>0; --i )
    {
        b_A->apply(y,tmp_result,Teuchos::NO_TRANS,1.0,0.0);
        y.update(1.0,tmp_result,d_coeffs[i-1],x,0.0);
    }

    // No iteration count, just return the polynomial order
    b_num_iters = num_coeffs-1;
}

} // namespace alea

