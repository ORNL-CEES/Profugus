//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SyntheticAcceleration.hh
 * \author Steven Hamilton
 * \brief  Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_MultiSplitting_hh
#define mc_solvers_MultiSplitting_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "LinearSystem_MultiSplitting.hh"
#include "LinearSystemFactory.hh"
#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class MultiSplitting
 * \brief Solver linear system using Multi-Splitting techniques
 */
//---------------------------------------------------------------------------//


class MultiSplitting 
{
  public:

    MultiSplitting(Teuchos::RCP<Teuchos::ParameterList>);
         
    Teuchos::RCP<const MATRIX> getMatrix() const { return d_A; }     

    // Implementation of solver
    void solve(Teuchos::RCP<MV> &x) const;

  private:

    Teuchos::RCP<const MATRIX> d_A;
    Teuchos::RCP<MV> d_b;

    // Parameter list
    Teuchos::RCP<Teuchos::ParameterList> d_pl;

    // Divergence tolerance
    MAGNITUDE d_divergence_tol;

    // Inner solver type
    std::string d_inner_solver;
    
    // Parameters for the generations of subdomains
    Teuchos::RCP<LinearSystem_MultiSplitting> d_multisplitting;
};

}

#endif // mc_solvers_SyntheticAcceleration_hh

