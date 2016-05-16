//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/AndersonSolver.t.hh
 * \author Steven Hamilton
 * \date   Wed Apr 01 11:01:28 2015
 * \brief  AndersonSolver template member definitions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_AndersonSolver_t_hh
#define SPn_solvers_AndersonSolver_t_hh

#include "comm/P_Stream.hh"
#include "AndersonSolver.hh"
#include "ThyraTraits.hh"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "NOX_StatusTest_RelativeNormF.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_Solver_AndersonAcceleration.H"
#include "NOX_Thyra_Group.H"

namespace profugus
{

//---------------------------------------------------------------------------//
/*\brief Constructor
 *
 * The following parameters can be specified on the ParameterList:
 *     tolerance(1e-4) - Stopping criterion on relative function norm
 *     max_itr(100)    - Number of iterations to perform before stopping
 * In addition, the NOX Anderson solver accepts the following parameters
 * on a nested "Anderson Parameters" list:
 *     "Storage Depth"                - Number of residual vectors retained
 *     "Mixing Parameter"             - Anderson mixing
 *     "Acceleration Start Iteration" - Perform specified Picard iterations
 *                                      before starting Anderson
 */
//---------------------------------------------------------------------------//
template <class T>
AndersonSolver<T>::AndersonSolver( Teuchos::RCP<const OP>               op,
                                   Teuchos::RCP<Teuchos::ParameterList> pl )
    : d_op(op)
    , d_pl(pl)
{
    REQUIRE(d_op != Teuchos::null);
    REQUIRE(d_pl != Teuchos::null);
}

//---------------------------------------------------------------------------//
//\brief Use NOX Anderson Acceleration to solve nonlinear problem
//---------------------------------------------------------------------------//
template <class T>
void AndersonSolver<T>::solve( Teuchos::RCP<MV> x )
{
    // Wrap operator into a Thyra ModelEvaluator
    d_model_eval = Teuchos::rcp( new ModelEvaluatorWrapper<T>(d_op) );
    CHECK( d_model_eval != Teuchos::null );

    // Create Thyra Vector from initial guess
    Teuchos::RCP<Thyra::MultiVectorBase<ST> > vec =
        ThyraTraits<T>::buildThyraMV(x,d_model_eval->get_x_space());
    CHECK( vec != Teuchos::null );

    // Build NOX "Group"
    // The ordinary constructor requires the model to support create_W_op,
    // so we pass in dummy arguments to force alternate constructor
    NOX::Thyra::Vector nox_vec(vec->col(0));
    Teuchos::RCP<NOX::Thyra::Group> nox_group(
        new NOX::Thyra::Group(nox_vec,d_model_eval,Teuchos::null,
            Teuchos::null,Teuchos::null,Teuchos::null) );
    CHECK( nox_group != Teuchos::null );

    // Status test
    // We probably want something more sophisticated long-term, but for
    // now just take a relative norm of the function with a maximum
    // iteration count.
    double tol = d_pl->get("tolerance",1e-4);
    Teuchos::RCP<NOX::StatusTest::Generic> status_rel_norm(
        new NOX::StatusTest::RelativeNormF(tol) );
    CHECK( status_rel_norm != Teuchos::null );

    int max_itr = d_pl->get("max_itr",100);
    Teuchos::RCP<NOX::StatusTest::Generic> status_max_itr(
        new NOX::StatusTest::MaxIters(max_itr) );
    CHECK( status_max_itr != Teuchos::null );

    Teuchos::RCP<NOX::StatusTest::Generic> status_combo(
        new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,status_rel_norm,
            status_max_itr) );
    CHECK( status_combo != Teuchos::null );

    // Set default mixing parameter to -1.0
    // This choice seems to almost always outperform +1.0
    Teuchos::sublist(d_pl,"Anderson Parameters")->get("Mixing Parameter",-1.0);

    // Build Anderson solver
    d_solver = Teuchos::rcp(
        new NOX::Solver::AndersonAcceleration(nox_group,status_combo,d_pl) );
    CHECK( d_solver != Teuchos::null );

    // Solve
    NOX::StatusTest::StatusType result = d_solver->solve();

    if( result == NOX::StatusTest::Converged )
        profugus::pout << "Anderson solver converged!" << profugus::endl;
    else
        profugus::pout << "Anderson solver did not converge!"
            << profugus::endl;

    // Copy solution back into x
    nox_vec.update(1.0,d_solver->getSolutionGroupPtr()->getX(),0.0);
}

} // end namespace profugus

#endif // SPn_solvers_AndersonSolver_t_hh

//---------------------------------------------------------------------------//
//                 end of AndersonSolver.t.hh
//---------------------------------------------------------------------------//
