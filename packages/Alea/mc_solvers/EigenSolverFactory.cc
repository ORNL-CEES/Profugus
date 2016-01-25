//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSolverFactory.cc
 * \author Massimiliano Lupo Pasini
 * \brief  Construct AleaSolver from ParameterList.
 */
//---------------------------------------------------------------------------//

#include <string>

#include "EigenSolverFactory.hh"

#include "AdditiveSchwarzWrapper.hh"
#include "PowerMethod.hh"
#include "MonteCarloEigenSolver.hh"
//#include "EigenMCSolver.hh"
//#include "EigenMCSA.hh"
//#include "EigenSequentialMC.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Public interface to build AleaSolver
 *
 * \param solver_type Name of linear solver to be constructed.
 * \param A Matrix for linear system.
 * \param pl ParameterList to be provided to solver.
 *
 *  Currently supported options for solver_type:
 *   - "power_method"  : PowerMethod
 *   - "monte_carlo"   : EigenMCSolver
 *   - "mcsa"          : EigenMCSA
 *   - "sequential_mc" : EigenSequentialMC
 *   - "none"        : No solver, this option is provided for when a solver is
 *                     being used as a preconditioner and no preconditioner is
 *                     requested.
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<AleaSolver>
EigenSolverFactory::buildSolver(std::string solver_type,
                                 Teuchos::RCP<const MATRIX> A,
                                 Teuchos::RCP<Teuchos::ParameterList> pl)
{
    VALIDATE(solver_type=="mcsa"                   ||
             solver_type=="sequential_mc"          ||
             solver_type=="monte_carlo"            ||
             solver_type=="power_method"           ||
             solver_type=="none",
            "Invalid solver_type.");

    // "mcsa" is shorthand for solver_type="synthetic_acceleration"
    //  with preconditioner="monte_carlo"
    /*if( solver_type == "mcsa" )
    {
	return Teuchos::rcp( new EigenMCSA(A,pl) );
    }

    if( solver_type == "sequential_mc" )
    {
        return Teuchos::rcp( new EigenSequentialMC(A,pl) );
    }
    */
    if( solver_type == "monte_carlo" )
    {
	Teuchos::RCP<MonteCarloEigenSolver> mc_solver( new MonteCarloEigenSolver(A,pl) );
        return Teuchos::rcp( new AdditiveSchwarzWrapper(A,mc_solver,pl) );
    }
    else if( solver_type == "power_method" )
    {
        return Teuchos::rcp( new PowerMethod(A,pl) );
    }

    // Return null if doesn't match any
    return Teuchos::null;
}

} // namespace alea

