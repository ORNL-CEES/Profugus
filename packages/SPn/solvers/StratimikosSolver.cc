//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/StratimikosSolver.cc
 * \author Chris Baker
 * \date   Thu May 17 14:18:09 2012
 * \brief  StratimikosSolver definitions
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include <Thyra_EpetraThyraWrappers.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "comm/global.hh"
#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "StratimikosSolver.hh"

#ifdef USE_MCLS
#include <MCLS_StratimikosAdapter.hpp>
#endif

namespace profugus
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
StratimikosSolver::StratimikosSolver(RCP_ParameterList db)
    : LinearSolver<MV,OP>(db)
{
    using Teuchos::sublist;

    // read database
    // get stratimikos database and convert to a ParameterList
    RCP_ParameterList builderplist = Teuchos::sublist(db, "Stratimikos");
    Check (!builderplist.is_null());

    if (db->isParameter("linear_solver_xml_file"))
    {
        Teuchos::updateParametersFromXmlFile(
            db->get<std::string>("linear_solver_xml_file"),
            builderplist.ptr());
    }

    if ( !Teuchos::isParameterType<std::string>(
             *builderplist, "Preconditioner Type") )
    {
        builderplist->set("Preconditioner Type", "None");
    }

    if ( !Teuchos::isParameterType<std::string>(
             *builderplist, "Linear Solver Type") )
    {
        builderplist->set("Linear Solver Type","AztecOO");
    }

    Stratimikos::DefaultLinearSolverBuilder builder;
#ifdef USE_MCLS
    MCLS::StratimikosAdapter<double>::setMCLSLinearSolveStrategyFactory(
	Teuchos::ptrFromRef(builder) );
    MCLS::StratimikosAdapter<double>::setMCLSPreconditioningStrategyFactory(
	Teuchos::ptrFromRef(builder) );
#endif
    builder.setParameterList(builderplist);
    b_label   = "Stratimikos " + builder.getLinearSolveStrategyName();
    d_factory = Thyra::createLinearSolveStrategy(builder);
    d_solver  = d_factory->createOp();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator and convert to Thyra operator
 *
 * \param A epetra operator
 */
void StratimikosSolver::set_operator(Teuchos::RCP<Epetra_Operator> A)
{
    // Create thyra operator
    d_thyraA = Thyra::epetraLinearOp(A);

    Ensure (d_thyraA != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set preconditioner.
 *
 * \param P epetra operator to be applied as preconditioner.
 */
void StratimikosSolver::set_preconditioner(Teuchos::RCP<Epetra_Operator> P)
{
    // Create thyra operator
    d_prec = Thyra::epetraLinearOp(P);

    Ensure( d_prec != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve linear system.
 *
 * \param ep_x pointer to solution vector
 * \param ep_b pointer to RHS vector
 */
void StratimikosSolver::solve(Teuchos::RCP<Epetra_MultiVector>       ep_x,
                              Teuchos::RCP<const Epetra_MultiVector> ep_b)
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ptrInArg;

    Insist (d_thyraA != Teuchos::null, "set_operator has not been called");
    Require (ep_x != Teuchos::null);
    Require (ep_b != Teuchos::null);
    Require (ep_x->NumVectors() == ep_b->NumVectors());
    Require (b_tolerance >= 0.0);
    Require (b_max_iters > 0);

    // Create Thyra preconditioner if we have an operator for it
    RCP_Prec prec;
    if( d_prec!=Teuchos::null )
    {
        prec = unspecifiedPrec(d_prec);
    }

    // Initialize LOWS
    if( prec != Teuchos::null )
    {
        Thyra::initializePreconditionedOp<double>(
            *d_factory, d_thyraA, prec, d_solver.ptr());
    }
    else
    {
        // Reuse as much information from previous solve as possible
        Thyra::initializeAndReuseOp<double>(
            *d_factory, d_thyraA, d_solver.ptr());
    }

    // make an Epetra view into the solution vector
    RCP<Thyra::MultiVectorBase<double> > x = Thyra::create_MultiVector(
        ep_x, d_thyraA->domain());

    // make an epetra vector view of the right-hand side
    RCP<const Thyra::MultiVectorBase<double> > b = Thyra::create_MultiVector(
        ep_b, d_thyraA->range());

    // solve
    Thyra::SolveCriteria<double> criteria;

    // measure convergence against |A*x-b| / |b|
    criteria.solveMeasureType = Thyra::SolveMeasureType(
        Thyra::SOLVE_MEASURE_NORM_RESIDUAL, Thyra::SOLVE_MEASURE_NORM_RHS);

    // set the tolerance
    criteria.requestedTol = b_tolerance;

    // set maximum iterations
    criteria.extraParameters = Teuchos::parameterList();
    criteria.extraParameters->set<int>("Maximum Iterations", b_max_iters);

    // solve
    Thyra::SolveStatus<double> solveStatus = Thyra::solve<double>(
        *d_solver, Thyra::NOTRANS, *b, x.ptr(), ptrInArg(criteria));

    // get the number of iterations for the solve
    b_num_iters = solveStatus.extraParameters->get<int>("Iteration Count");

    double final_res = solveStatus.achievedTol;

    Thyra::ESolveStatus status = solveStatus.solveStatus;

    if( b_verbosity >=
            LinearSolver<Epetra_MultiVector,Epetra_Operator>::LOW )
    {
        profugus::pout << b_label << " performed " << b_num_iters
            << " iterations.  Final residual norm is "
            << profugus::scientific << profugus::setprecision(3)
            << final_res << profugus::endl;
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of StratimikosSolver.cc
//---------------------------------------------------------------------------//
