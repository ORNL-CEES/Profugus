//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/LinearSolverFactory.hh
 * \author Steven Hamilton
 * \brief  Construct Tpetra_CrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_LinearSolverFactory_hh
#define Alea_mc_solvers_LinearSolverFactory_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"
#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class LinearSolverFactory
 * \brief Construct linear solver from MATRIX and ParameterList.
 */
//---------------------------------------------------------------------------//
class LinearSolverFactory
{
  private:

    // Pure static, disallow construction
    LinearSolverFactory(){};

  public:

    static Teuchos::RCP<AleaSolver> buildSolver(
        std::string solver_type,
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl );

};

}

#endif // Alea_mc_solvers_LinearSolverFactory_hh

