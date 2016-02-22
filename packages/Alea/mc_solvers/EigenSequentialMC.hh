//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSequentialMC.hh
 * \author Massimiliano Lupo Pasini
 * \brief  Perform Sequential Monte Carlo for dominant eigenvalue.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_EigenSequentialMC_hh
#define mc_solvers_EigenSequentialMC_hh

#include "AleaSolver.hh"
#include "MonteCarloEigenSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class EigenSequentialMC
 * \brief Solve eigen problem using power iteration.
 *
 * This class solves an eigenvalue problem of equations using power iterations.
 */
//---------------------------------------------------------------------------//
class EigenSequentialMC : public AleaSolver
{
  public:

    EigenSequentialMC(Teuchos::RCP<const MATRIX> A,
                        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    SCALAR d_divergence_tol;

    Teuchos::RCP<MonteCarloEigenSolver> d_mc_solver;
};

}

#endif // mc_solvers_EigenSequentialMC_hh

