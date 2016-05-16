//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/RichardsonIteration.hh
 * \author Steven Hamilton
 * \brief  Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_RichardsonIteration_hh
#define Alea_mc_solvers_RichardsonIteration_hh

#include "AleaSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class RichardsonIteration
 * \brief Solve linear system using Richardson iteration.
 *
 * This class solves a linear system of equations using Richardson iteration
 * with an optional preconditioner.
 * Richardson iteration for solving the system
 * \f$ \textbf{Ax} = \textbf{b} \f$ is given by:
 *
 * \f{eqnarray*}{
 *     \textbf{r}^{k} &=& \textbf{b} - \textbf{Ax}^{k} \\
 *     \textbf{x}^{k+1} &=& \textbf{x}^{k} +
 *        \textbf{M}^{-1}\textbf{r}^{k},
 * \f}
 * where \f$\textbf{M}^{-1}\f$ is a preconditioner.
 */
//---------------------------------------------------------------------------//
class RichardsonIteration : public AleaSolver
{
  public:

    RichardsonIteration(Teuchos::RCP<const MATRIX> A,
                        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    Teuchos::RCP<const AleaSolver> d_P;

    SCALAR d_divergence_tol;

};

}

#endif // Alea_mc_solvers_RichardsonIteration_hh

