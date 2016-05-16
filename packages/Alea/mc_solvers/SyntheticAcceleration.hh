//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/SyntheticAcceleration.hh
 * \author Steven Hamilton
 * \brief  Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_SyntheticAcceleration_hh
#define Alea_mc_solvers_SyntheticAcceleration_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class SyntheticAcceleration
 * \brief Solver linear system using synthetic acceleration
 *
 * The synthetic acceleration algorithm for solving the system
 * \f$ \textbf{Ax} = \textbf{b} \f$ is given by:
 *
 * \f{eqnarray*}{
 *     \textbf{r}^{k} &=& \textbf{b} - \textbf{Ax}^{k} \\
 *     \textbf{x}^{k+1/2} &=& \textbf{x}^{k} + \textbf{r}^{k} \\
 *     \textbf{r}^{k+1/2} &=& \textbf{b} - \textbf{Ax}^{k+1/2} \\
 *     \textbf{x}^{k+1} &=& \textbf{x}^{k+1/2} +
 *        \textbf{M}^{-1}\textbf{r}^{k+1/2},
 * \f}
 * where \f$\textbf{M}^{-1}\f$ indicates an approximate solution of a linear
 * system (or solution of an approximate linear system) using some
 * preconditioner  In this class, the preconditioner can be provided by
 * any AleaSolver class (see LinearSolverFactory for details).
 * In particular, if the MonteCarloSolver class is used as a preconditioner,
 * then the Monte Carlo synthetic acceleration (MCSA) algorithm results.
 */
//---------------------------------------------------------------------------//
class SyntheticAcceleration : public AleaSolver
{
  public:

    SyntheticAcceleration(Teuchos::RCP<const MATRIX> A,
         Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    // Preconditioner
    Teuchos::RCP<AleaSolver> d_preconditioner;

    // Damping factor
    SCALAR d_damping;

    // Divergence tolerance
    MAGNITUDE d_divergence_tol;
};

}

#endif // Alea_mc_solvers_SyntheticAcceleration_hh

