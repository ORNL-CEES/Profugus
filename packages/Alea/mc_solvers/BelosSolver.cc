//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BelosSolver.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "BelosSolver.hh"
#include "LinearSolverFactory.hh"

#include "BelosSolverFactory.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by following PL entries on the nested "Belos"
 * sublist:
 *  - belos_type(string)                : see Belos documentation for valid
 *                                        options ("Block GMRES")
 *  - preconditioner(string)            : any MCREX solver ("none")
 *  - Maximum Iterations(int)           : >0 (1000)
 *  - Convergence Tolerance (MAGNITUDE) : >0.0 (1.0e-6)
 *  - verbosity(string)                 : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
BelosSolver::BelosSolver( Teuchos::RCP<const MATRIX> A,
                          Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get belos pl
    Teuchos::RCP<Teuchos::ParameterList> belos_pl =
        Teuchos::sublist(pl,"Belos");

    // Override default parameters if present on sublist
    this->setParameters(belos_pl);

    // Set default values on PL
    if( !belos_pl->isType<MAGNITUDE>("Convergence Tolerance") )
    {
        MAGNITUDE tol = 1.0e-6;
        belos_pl->set("Convergence Tolerance",tol);
    }

    if( !belos_pl->isType<LO>("Maximum Iterations") )
    {
        LO max_itr = 1000;
        belos_pl->set("Maximum Iterations",max_itr);
    }

    // Build preconditioner, for now we just use LinearSolverFactory
    std::string prec_type = pl->get("preconditioner","none");
    d_prec = LinearSolverFactory::buildSolver(prec_type,b_A,pl);
    if( prec_type.compare("monte_carlo") == 0 ||
        prec_type.compare("mcsa") == 0 )
    {
        belos_pl->set("Flexible Gmres",true);
    }

    // Set verbosity -- always enable warnings
    LO verb = Belos::Warnings;
    if( b_verbosity >= LOW )
    {
        verb += Belos::FinalSummary;
    }
    if( b_verbosity >= HIGH )
    {
        verb += Belos::StatusTestDetails;
    }
    belos_pl->set("Verbosity",verb);

    // Build Belos solver
    std::string belos_type = belos_pl->get("belos_type","Block GMRES");
    Belos::SolverFactory<SCALAR,MV,OP> factory;
    d_belos_solver = factory.create(belos_type,belos_pl);

    b_label = "BelosSolver";
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Solve linear system using Belos.
 */
//---------------------------------------------------------------------------//
void BelosSolver::applyImpl(const MV &x, MV &y) const
{
    // Create linear problem
    Teuchos::RCP<Belos::LinearProblem<SCALAR,MV,OP> > problem(
            new Belos::LinearProblem<SCALAR,MV,OP>() );
    problem->setOperator(b_A);
    problem->setLHS(Teuchos::rcpFromRef(y));
    problem->setRHS(Teuchos::rcpFromRef(x));
    if( d_prec != Teuchos::null )
        problem->setRightPrec(d_prec);
    problem->setProblem();

    d_belos_solver->setProblem(problem);

    // Solve
    Belos::ReturnType result = d_belos_solver->solve();

    b_num_iters = d_belos_solver->getNumIters();
    if( b_verbosity >= LOW )
    {
        std::cout << b_label
                  << (result == Belos::Converged ? " converged "
                                                 : " did not converge ")
                  << "after " << b_num_iters << " iterations." << std::endl;
    }
}

} // namespace alea

