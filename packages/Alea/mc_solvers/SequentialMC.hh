//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SequentialMC.hh
 * \author Steven Hamilton
 * \brief  Performs Sequential Monte Carlo.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_SequentialMC_hh
#define mc_solvers_SequentialMC_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class SequentialMC
 * \brief Solver linear system using Sequential Monte Carlo
 *
 * The Sequential Monte Carlo algorithm for solving the system
 * \f$ \textbf{Ax} = \textbf{b} \f$ is given by:
 *
 * \f{eqnarray*}{
 *     \textbf{r}^{k} &=& \textbf{b} - \textbf{Ax}^{k} \\
 *    \textbf{dx} &=&  \textbf{M}^{-1}\textbf{r}^{k} \\
 *    \textbf{x}^{k+1} &=& \textbf{x}^{k} + \textbf{dx},
 * \f}
 * where \f$\textbf{M}^{-1}\f$ indicates an approximate solution of a linear
 * system (or solution of an approximate linear system) using some
 * preconditioner  In this class, the preconditioner is provided by 
 * the MonteCarloSolver class,
 * giving birth to the so called Sequential Monte Carlo.
 */
//---------------------------------------------------------------------------//
class SequentialMC : public AleaSolver
{
  public:

    SequentialMC(Teuchos::RCP<const MATRIX> A,
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

#endif // mc_solvers_SyntheticAcceleration_hh

