//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSolverFactory.hh
 * \author Massimiliano Lupo Pasini
 * \brief  Construct Tpetra_CrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_EigenSolverFactory_hh
#define mc_solvers_EigenSolverFactory_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"
#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class EigenSolverFactory
 * \brief Construct linear solver from MATRIX and ParameterList.
 */
//---------------------------------------------------------------------------//
class EigenSolverFactory
{
  private:

    // Pure static, disallow construction
    EigenSolverFactory(){};

  public:

    static Teuchos::RCP<AleaSolver> buildSolver(
        std::string solver_type,
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl );

};

}

#endif // mc_solvers_EigenSolverFactory_hh

