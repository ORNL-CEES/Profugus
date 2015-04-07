//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/StratimikosSolver.t.hh
 * \author Chris Baker
 * \date   Thu May 17 14:18:09 2012
 * \brief  StratimikosSolver template member definitions
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_StratimikosSolver_t_hh
#define solvers_StratimikosSolver_t_hh

#include <iostream>

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "comm/global.hh"
#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "StratimikosSolver.hh"
#include "LinAlgTypedefs.hh"
#include "ThyraTraits.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
template <class T>
StratimikosSolver<T>::StratimikosSolver(Teuchos::RCP<ParameterList> db)
    : LinearSolver<T>(db)
    , d_updated_operator( false )
{
    using Teuchos::sublist;

    // read database
    // get stratimikos database and convert to a ParameterList
    Teuchos::RCP<ParameterList> builderplist =
        Teuchos::sublist(db, "Stratimikos");
    CHECK(!builderplist.is_null());

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
    builder.setParameterList(builderplist);
    b_label   = "Stratimikos " + builder.getLinearSolveStrategyName();
    d_factory = Thyra::createLinearSolveStrategy(builder);
    d_solver  = d_factory->createOp();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
//
//---------------------------------------------------------------------------//
/*!
 * \brief Solve linear system.
 *
 * \param x pointer to solution vector
 * \param b pointer to RHS vector
 */
template <class T>
void StratimikosSolver<T>::solve(Teuchos::RCP<MV>       x,
                                 Teuchos::RCP<const MV> b)
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ptrInArg;

    INSIST(d_thyraA != Teuchos::null, "set_operator has not been called");
    REQUIRE(x != Teuchos::null);
    REQUIRE(b != Teuchos::null);
    REQUIRE(MVT::GetNumberVecs(*x) == MVT::GetNumberVecs(*b));
    REQUIRE(b_tolerance >= 0.0);
    REQUIRE(b_max_iters > 0);

    // Create Thyra preconditioner if we have an operator for it
    RCP<const Prec> prec;
    if( d_prec!=Teuchos::null )
    {
        prec = unspecifiedPrec(d_prec);
    }

    // Initialize LOWS
    if( prec != Teuchos::null )
    {
        Thyra::initializePreconditionedOp<ST>(
            *d_factory, d_thyraA, prec, d_solver.ptr());
    }
    // If the operator has changed but we are reusing the preconditioner then
    // reinitialize the linearOpWithSolve object.
    else if ( d_updated_operator )
    {
        // Reuse as much information from previous solve as possible
        Thyra::initializeAndReuseOp<ST>(
            *d_factory, d_thyraA, d_solver.ptr());
        d_updated_operator = false;
    }

    // make Thyra view into the solution vector
    RCP<Thyra::MultiVectorBase<ST> > thyra_x =
        ThyraTraits<T>::buildThyraMV(x,d_thyraA->domain());

    // make Thyra view of the right-hand side
    RCP<const Thyra::MultiVectorBase<ST> > thyra_b =
        ThyraTraits<T>::buildThyraConstMV(b,d_thyraA->range());

    // solve
    Thyra::SolveCriteria<ST> criteria;

    // measure convergence against |A*x-b| / |b|
    criteria.solveMeasureType = Thyra::SolveMeasureType(
        Thyra::SOLVE_MEASURE_NORM_RESIDUAL, Thyra::SOLVE_MEASURE_NORM_RHS);

    // set the tolerance
    criteria.requestedTol = b_tolerance;

    // set maximum iterations
    criteria.extraParameters = Teuchos::parameterList();
    criteria.extraParameters->template set<int>("Maximum Iterations", b_max_iters);

    // solve
    Thyra::SolveStatus<ST> solveStatus = Thyra::solve<ST>(
        *d_solver, Thyra::NOTRANS, *thyra_b, thyra_x.ptr(), ptrInArg(criteria));

    // get the number of iterations for the solve
    b_num_iters = solveStatus.extraParameters->template get<int>("Iteration Count");

    if( b_verbosity >= LinearSolver<T>::LOW )
    {
        ST final_res = solveStatus.achievedTol;

        profugus::pout << b_label << " performed " << b_num_iters
            << " iterations.  Final residual norm is "
            << profugus::scientific << profugus::setprecision(3)
            << final_res << profugus::endl;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set operator and convert to Thyra operator
 *
 * \param A operator
 */
//---------------------------------------------------------------------------//
template <class T>
void StratimikosSolver<T>::set_operator(Teuchos::RCP<const OP> A)
{
    // Create thyra operator
    d_thyraA = ThyraTraits<T>::buildThyraOP(A);

    // Indicate that the operator has been updated.
    d_updated_operator = true;

    ENSURE(d_thyraA != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set preconditioner.
 *
 * \param P operator to be applied as preconditioner.
 */
//---------------------------------------------------------------------------//
template <class T>
void StratimikosSolver<T>::set_preconditioner(Teuchos::RCP<const OP> P)
{
    // Create thyra operator
    d_prec = ThyraTraits<T>::buildThyraOP(P);

    ENSURE( d_prec != Teuchos::null );
}

} // end namespace profugus

#endif //solvers_StratimikosSolver_t_hh

//---------------------------------------------------------------------------//
//                 end of StratimikosSolver.t.hh
//---------------------------------------------------------------------------//
