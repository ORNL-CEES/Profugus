//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BelosSolver.cc
 * \author Steven Hamilton
 * \brief  BelosSolver template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_BelosSolver_t_hh
#define solvers_BelosSolver_t_hh

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
    }
    belos_pl->set("Verbosity",verb);

    // Build Belos solver
    std::string belos_type = belos_pl->get("belos_type","Block GMRES");
    Belos::SolverFactory<ST,MV,OP> factory;
    d_belos_solver = factory.create(belos_type,belos_pl);

    b_label = "Belos" + belos_type;
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
    if( b_verbosity >= LinearSolver<T>::LOW )
    {
        std::cout << b_label
                  << (result == Belos::Converged ? " converged "
                                                 : " did not converge ")
                  << "after " << b_num_iters << " iterations." << std::endl;
    }
}

} // namespace profugus

#endif // solvers_BelosSolver_t_hh

