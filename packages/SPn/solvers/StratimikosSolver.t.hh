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
#include "TpetraTypedefs.hh"

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
template <class MV, class OP>
StratimikosSolver<MV,OP>::StratimikosSolver(RCP_ParameterList db)
    : LinearSolver<MV,OP>(db)
    , d_updated_operator( false )
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
//
//---------------------------------------------------------------------------//
/*!
 * \brief Solve linear system.
 *
 * \param x pointer to solution vector
 * \param b pointer to RHS vector
 */
template <class MV, class OP>
void StratimikosSolver<MV,OP>::solve(Teuchos::RCP<MV>       x,
                                     Teuchos::RCP<const MV> b)
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ptrInArg;

    Insist (d_thyraA != Teuchos::null, "set_operator has not been called");
    Require (x != Teuchos::null);
    Require (b != Teuchos::null);
    Require (MVT::GetNumberVecs(*x) == MVT::GetNumberVecs(*b));
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
    // If the operator has changed but we are reusing the preconditioner then
    // reinitialize the linearOpWithSolve object.
    else if ( d_updated_operator )
    {
        // Reuse as much information from previous solve as possible
        Thyra::initializeAndReuseOp<double>(
            *d_factory, d_thyraA, d_solver.ptr());
        d_updated_operator = false;
    }

    // make Thyra view into the solution vector
    RCP<Thyra::MultiVectorBase<double> > thyra_x =
        buildThyraMV(x,d_thyraA->domain());

    // make Thyra view of the right-hand side
    RCP<const Thyra::MultiVectorBase<double> > thyra_b =
        buildThyraConstMV(b,d_thyraA->range());

    // solve
    Thyra::SolveCriteria<double> criteria;

    std::vector<double> b_norm(1);
    MVT::MvNorm(*b,b_norm);

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
        *d_solver, Thyra::NOTRANS, *thyra_b, thyra_x.ptr(), ptrInArg(criteria));

    std::vector<double> x_norm(1);
    MVT::MvNorm(*x,x_norm);

    // get the number of iterations for the solve
    b_num_iters = solveStatus.extraParameters->get<int>("Iteration Count");

    double final_res = solveStatus.achievedTol;

    Thyra::ESolveStatus status = solveStatus.solveStatus;

    if( b_verbosity >= LinearSolver<MV,OP>::LOW )
    {
        profugus::pout << b_label << " performed " << b_num_iters
            << " iterations.  Final residual norm is "
            << profugus::scientific << profugus::setprecision(3)
            << final_res << profugus::endl;
    }
}

//---------------------------------------------------------------------------//
// TEMPLATE SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set operator and convert to Thyra operator
 *
 * \param A operator
 */
template <>
void StratimikosSolver<Epetra_MultiVector,Epetra_Operator>::set_operator(
    Teuchos::RCP<Epetra_Operator> A)
{
    // Create thyra operator
    d_thyraA = Thyra::epetraLinearOp(A);

    // Indicate that the operator has been updated.
    d_updated_operator = true;

    Ensure (d_thyraA != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set preconditioner.
 *
 * \param P operator to be applied as preconditioner.
 */
template <>
void StratimikosSolver<Epetra_MultiVector,Epetra_Operator>::set_preconditioner(
    Teuchos::RCP<Epetra_Operator> P)
{
    // Create thyra operator
    d_prec = Thyra::epetraLinearOp(P);

    Ensure( d_prec != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wrap MV into a Thyra vector
 */
template <>
Teuchos::RCP<Thyra::MultiVectorBase<double> >
StratimikosSolver<Epetra_MultiVector,Epetra_Operator>::buildThyraMV(
    Teuchos::RCP<Epetra_MultiVector> x,
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > space) const
{
    return Thyra::create_MultiVector(x,space);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wrap const MV into a Thyra vector
 */
template <>
Teuchos::RCP<const Thyra::MultiVectorBase<double> >
StratimikosSolver<Epetra_MultiVector,Epetra_Operator>::buildThyraConstMV(
    Teuchos::RCP<const Epetra_MultiVector> x,
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > space) const
{
    return Thyra::create_MultiVector(x,space);
}

/*!
 * \brief Set operator and convert to Thyra operator
 *
 * \param A operator
 */
template <>
void StratimikosSolver<
    Tpetra::MultiVector<double,int,int,KokkosClassic::SerialNode>,
    Tpetra::Operator<double,int,int,KokkosClassic::SerialNode> >::set_operator(
        Teuchos::RCP<Tpetra::Operator<double,int,int,KokkosClassic::SerialNode> > A)
{
    // Create thyra operator
    Teuchos::RCP<Thyra::VectorSpaceBase<double> > rangeSpace =
        Thyra::tpetraVectorSpace<double>(A->getRangeMap());
    Teuchos::RCP<Thyra::VectorSpaceBase<double> > domainSpace =
        Thyra::tpetraVectorSpace<double>(A->getDomainMap());
    d_thyraA = Thyra::tpetraLinearOp<double,int,int,KokkosClassic::SerialNode>(
        rangeSpace,domainSpace,A);

    // Indicate that the operator has been updated.
    d_updated_operator = true;

    Ensure (d_thyraA != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set preconditioner.
 *
 * \param P operator to be applied as preconditioner.
 */
template <>
void
StratimikosSolver<Tpetra_MultiVector,Tpetra_Operator>::set_preconditioner(
    Teuchos::RCP<Tpetra_Operator> P)
{
    // Create thyra operator
    Teuchos::RCP<const Thyra::VectorSpaceBase<SCALAR> > rangeSpace =
        Thyra::tpetraVectorSpace<SCALAR>(P->getRangeMap());
    Teuchos::RCP<const Thyra::VectorSpaceBase<SCALAR> > domainSpace =
        Thyra::tpetraVectorSpace<SCALAR>(P->getDomainMap());
    d_prec = Thyra::tpetraLinearOp<SCALAR,LO,GO,NODE>(
        rangeSpace,domainSpace,P);

    Ensure( d_prec != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wrap MV into a Thyra vector
 */
template <>
Teuchos::RCP<Thyra::MultiVectorBase<SCALAR> >
StratimikosSolver<Tpetra_MultiVector,Tpetra_Operator>::buildThyraMV(
    Teuchos::RCP<Tpetra_MultiVector> x,
    Teuchos::RCP<const Thyra::VectorSpaceBase<SCALAR> > space) const
{
    return Thyra::createMultiVector(x,space);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wrap const MV into a Thyra vector
 */
template <>
Teuchos::RCP<const Thyra::MultiVectorBase<SCALAR> >
StratimikosSolver<Tpetra_MultiVector,Tpetra_Operator>::buildThyraConstMV(
    Teuchos::RCP<const Tpetra_MultiVector> x,
    Teuchos::RCP<const Thyra::VectorSpaceBase<SCALAR> > space) const
{
    return Thyra::createConstMultiVector(x,space);
}

} // end namespace profugus

#endif //solvers_StratimikosSolver_t_hh

//---------------------------------------------------------------------------//
//                 end of StratimikosSolver.t.hh
//---------------------------------------------------------------------------//
