//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolverFactory.cc
 * \author Steven Hamilton
 * \brief  Construct AleaSolver from ParameterList.
 */
//---------------------------------------------------------------------------//

#include <string>

#include "LinearSolverFactory.hh"

#include "AdditiveSchwarzWrapper.hh"
#include "BelosSolver.hh"
#include "ChebyshevIteration.hh"
//#include "MonteCarloSolver.hh"
#include "PolynomialPreconditioner.hh"
#include "RichardsonIteration.hh"
#include "SyntheticAcceleration.hh"

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
 *   - "richardson"  : RichardsonIteration
 *   - "monte_carlo" : MonteCarloSolver
 *   - "mcsa"        : SyntheticAcceleration with MonteCarloSolver
 *                     preconditioner.
 *   - "belos"       : BelosSolver
 *   - "chebyshev"   : ChebyshevIteration
 *   - "synthetic_acceleration" : SyntheticAcceleration
 *   - "polynomial"  : PolynomialPreconditioner (only useful as preconditioner)
 *   - "none"        : No solver, this option is provided for when a solver is
 *                     being used as a preconditioner and no preconditioner is
 *                     requested.
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<AleaSolver>
LinearSolverFactory::buildSolver(std::string solver_type,
                                 Teuchos::RCP<const MATRIX> A,
                                 Teuchos::RCP<Teuchos::ParameterList> pl)
{
    TEUCHOS_ASSERT( solver_type=="belos"                  ||
                    solver_type=="chebyshev"              ||
                    solver_type=="mcsa"                   ||
                    solver_type=="monte_carlo"            ||
                    solver_type=="polynomial"             ||
                    solver_type=="richardson"             ||
                    solver_type=="synthetic_acceleration" ||
                    solver_type=="none" );

    // "mcsa" is shorthand for solver_type="synthetic_acceleration"
    //  with preconditioner="monte_carlo"
    if( solver_type == "mcsa" )
    {
        solver_type = "synthetic_acceleration";
        pl->set("preconditioner","monte_carlo");
    }

    if( solver_type == "richardson" )
    {
        return Teuchos::rcp( new RichardsonIteration(A,pl) );
    }
    else if( solver_type == "monte_carlo" )
    {
        /*
        // The Monte Carlo solver is a local domain solver only
        // We wrap it into an AdditiveSchwarz to give a global solver
        Teuchos::RCP<MonteCarloSolver> mc_solver( new MonteCarloSolver(A,pl) );

        return Teuchos::rcp( new AdditiveSchwarzWrapper(A,mc_solver,pl) );
        */
        TEUCHOS_ASSERT(false);
    }
    else if( solver_type == "synthetic_acceleration" )
    {
        return Teuchos::rcp( new SyntheticAcceleration(A,pl) );
    }
    else if( solver_type == "belos" )
    {
        return Teuchos::rcp( new BelosSolver(A,pl) );
    }
    else if( solver_type == "chebyshev" )
    {
        return Teuchos::rcp( new ChebyshevIteration(A,pl) );
    }
    else if( solver_type == "polynomial" )
    {
        return Teuchos::rcp( new PolynomialPreconditioner(A,pl) );
    }

    // Return null if doesn't match any
    return Teuchos::null;
}

} // namespace alea

