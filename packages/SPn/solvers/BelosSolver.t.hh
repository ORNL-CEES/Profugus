//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/BelosSolver.t.hh
 * \author Steven Hamilton
 * \brief  BelosSolver template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_BelosSolver_t_hh
#define SPn_solvers_BelosSolver_t_hh

#include <iterator>
#include <string>

#include "BelosSolver.hh"

#include "BelosSolverFactory.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace profugus
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
 *  - Maximum Iterations(int)           : >0 (1000)
 *  - Convergence Tolerance (MAGNITUDE) : >0.0 (1.0e-6)
 *  - verbosity(string)                 : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
template <class T>
BelosSolver<T>::BelosSolver( Teuchos::RCP<Teuchos::ParameterList> pl )
  : LinearSolver<T>(pl)
{
    // Get belos pl
    Teuchos::RCP<Teuchos::ParameterList> belos_pl =
        Teuchos::sublist(pl,"Belos");

    // Set stopping criteria on Belos pl
    if( !belos_pl->isType<ST>("Convergence Tolerance") )
    {
        belos_pl->set("Convergence Tolerance",b_tolerance);
    }

    if( !belos_pl->isType<LO>("Maximum Iterations") )
    {
        belos_pl->set("Maximum Iterations",b_max_iters);
    }

    // Set verbosity -- always enable warnings
    LO verb = Belos::Warnings;
    if( b_verbosity >= LinearSolver<T>::LOW )
    {
        verb += Belos::FinalSummary;
    }
    if( b_verbosity >= LinearSolver<T>::HIGH )
    {
        verb += Belos::StatusTestDetails;
        belos_pl->set("Output Frequency",1);
    }
    belos_pl->set("Verbosity",verb);

    d_belos_type = b_db->get("belos_type","Pseudo Block GMRES");
    d_belos_pl = belos_pl;

    // Set residual scaling based on RHS rather than initial
    // residual, which is the default behavior
    d_belos_pl->set("Implicit Residual Scaling","Norm of RHS");
    d_belos_pl->set("Explicit Residual Scaling","Norm of RHS");

    // Build Belos solver
    Belos::SolverFactory<ST,MV,OP> factory;
    d_belos_solver = factory.create(d_belos_type,d_belos_pl);

    b_label = "Belos " + d_belos_type;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve linear system using Belos.
 */
//---------------------------------------------------------------------------//
template <class T>
void BelosSolver<T>::solve(Teuchos::RCP<MV> x, Teuchos::RCP<const MV> b)
{
    INSIST( b_A != Teuchos::null, "set_operator has not been called");

    // Make sure solver has latest parameters
    d_belos_solver->setParameters(d_belos_pl);

    // Create linear problem
    Teuchos::RCP<Belos::LinearProblem<ST,MV,OP> > problem(
            new Belos::LinearProblem<ST,MV,OP>() );
    problem->setOperator(b_A);
    problem->setLHS(x);
    problem->setRHS(b);
    if( d_P != Teuchos::null )
        problem->setRightPrec(d_P);
    problem->setProblem();

    d_belos_solver->setProblem(problem);

    // Solve
    Belos::ReturnType result = d_belos_solver->solve();

    b_num_iters = d_belos_solver->getNumIters();
    b_converged = (result == Belos::Converged);
    if( b_verbosity >= LinearSolver<T>::LOW )
    {
        std::cout << b_label
                  << (b_converged ? " converged "
                                  : " did not converge ")
                  << "after " << b_num_iters << " iterations." << std::endl;
    }
}

} // namespace profugus

#endif // SPn_solvers_BelosSolver_t_hh

