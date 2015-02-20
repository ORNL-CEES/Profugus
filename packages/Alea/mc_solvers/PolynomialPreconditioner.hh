//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialPreconditioner.hh
 * \author Steven Hamilton
 * \brief  Perform Chebyshev iteration.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_PolynomialPreconditioner_hh
#define mc_solvers_PolynomialPreconditioner_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class PolynomialPreconditioner
 * \brief Apply matrix polynomial.
 *
 * This class generates a matrix polynomial (using a Polynomial
 * constructed by the PolynomialFactory class)
 * and applies it to a vector using Horner's method.  This class is intended
 * to be used as a preconditioner to a linear solver.
 */
//---------------------------------------------------------------------------//
class PolynomialPreconditioner : public AleaSolver
{
  public:

    PolynomialPreconditioner(Teuchos::RCP<const MATRIX> A,
                             Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    Teuchos::ArrayRCP<const SCALAR> d_coeffs;
};

}

#endif // mc_solvers_PolynomialPreconditioner_hh

