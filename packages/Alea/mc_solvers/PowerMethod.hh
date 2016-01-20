//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerMethod.hh
 * \author Massimiliano Lupo Pasini
 * \brief  Perform power iteration.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_PowerMethod_hh
#define mc_solvers_PowerMethod_hh

#include "AleaSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class PowerMethod
 * \brief Solve eigen problem using power iteration.
 *
 * This class solves an eigenvalue problem of equations using power iterations.
 */
//---------------------------------------------------------------------------//
class PowerMethod : public AleaSolver
{
  public:

    PowerMethod(Teuchos::RCP<const MATRIX> A,
                        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    SCALAR d_divergence_tol;

};

}

#endif // mc_solvers_PowerMethod_hh

