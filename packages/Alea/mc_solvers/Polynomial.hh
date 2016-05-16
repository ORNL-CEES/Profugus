//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/Polynomial.hh
 * \author Steven Hamilton
 * \brief  Polynomial class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_Polynomial_hh
#define Alea_mc_solvers_Polynomial_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AnasaziTypes.hpp"

#include "AleaTypedefs.hh"

#include "PolynomialBasis.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class Polynomial
 * \brief Base class for polynomial types.
 */
//---------------------------------------------------------------------------//
class Polynomial
{
  public:

    // Access coefficient values.
    virtual Teuchos::ArrayRCP<const SCALAR>
    getCoeffs(const PolynomialBasis &target = PolynomialBasis("power")) const;

    //! \brief Access PolynomialBasis object used internally.
    virtual Teuchos::RCP<const PolynomialBasis> getNativeBasis() const
    {
        REQUIRE( !b_native_basis.is_null() );
        return b_native_basis;
    }

    virtual ~Polynomial(){};

  protected:

    // Constructor, protected to avoid direct construction
    Polynomial(Teuchos::RCP<const MATRIX> A,
               Teuchos::RCP<Teuchos::ParameterList> pl);

    enum Verbosity_Level { NONE, LOW, MEDIUM, HIGH };

    // Original matrix
    Teuchos::RCP<const MATRIX> b_A;

    // ParameterLists
    Teuchos::RCP<Teuchos::ParameterList> b_pl;
    Teuchos::RCP<Teuchos::ParameterList> b_poly_pl;

    // Polynomial order
    int b_m;

    // Polynomial coefficients
    Teuchos::ArrayRCP<SCALAR> b_coeffs;

    // Polynomial basis
    Teuchos::RCP<PolynomialBasis> b_native_basis;

    // Output level
    Verbosity_Level b_verbosity;
};

}

#endif // Alea_mc_solvers_Polynomial_hh

