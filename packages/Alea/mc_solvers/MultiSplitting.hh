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

  private:

    Teuchos::RCP<MATRIX> d_A;
    Teuchos::RCP<MV> d_b;

    // Implementation of apply
    void solve(const MV &x, MV &y) const;

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

