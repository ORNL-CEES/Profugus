//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/Polynomial.cc
 * \author Steven Hamilton
 * \brief  Polynomial class definitions.
 */
//---------------------------------------------------------------------------//

#include "comm/global.hh"

#include "Polynomial.hh"
#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
//---------------------------------------------------------------------------//
Polynomial::Polynomial(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<Teuchos::ParameterList> pl)
  : b_A(A)
  , b_pl(pl)
{
    REQUIRE( b_pl != Teuchos::null );

    b_poly_pl = Teuchos::sublist(b_pl,"Polynomial");

    b_m = b_poly_pl->get<int>("polynomial_order",1);

    // Set verbosity for polynomial construction
    std::string verb = b_poly_pl->get("verbosity","none");
    if( verb == "high" )
        b_verbosity = HIGH;
    else if( verb == "medium" )
        b_verbosity = MEDIUM;
    else if( verb == "low" )
        b_verbosity = LOW;
    else
        b_verbosity = NONE;

    // Silence processors other than 0
    if( profugus::node() != 0 )
        b_verbosity = NONE;
}

//---------------------------------------------------------------------------//
// PUBLIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Access polynomial coefficients.
 *
 * Each polynomial type constructs a set of coefficients \f$ c_i \f$ such
 * that \f$ p(x) = \sum_{i=0}^m c_i b_\mathrm{native}(x)^i \f$, where
 * \f$ b_\mathrm{native}(x) = \alpha + \beta x \f$ is a linear basis that is
 * chosen by the particular polynomial type.  The calling function
 * (or the user) will typically want to work in a prescribed ``target''
 * basis.  This function returns the coefficients with respect to that
 * target basis.
 */
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<const SCALAR>
Polynomial::getCoeffs(const PolynomialBasis &target) const
{
    REQUIRE( b_native_basis != Teuchos::null );
    REQUIRE( b_coeffs != Teuchos::null );
    REQUIRE( b_coeffs.size() == (b_m+1) );

    // Convert coefficients from native basis to target
    Teuchos::ArrayRCP<SCALAR> target_coeffs;
    target_coeffs = target.transformBasis(b_coeffs,*b_native_basis);
    REQUIRE( target_coeffs.size() == (b_m+1) );

    if( b_verbosity >= MEDIUM )
    {
        std::cout << "Polynomial Coefficients: " << std::endl;
        for( int i=0; i<b_m+1; ++i )
        {
            std::cout << i << " " << b_coeffs[i] << std::endl;
        }
    }

    return target_coeffs;
}

} // namespace alea

