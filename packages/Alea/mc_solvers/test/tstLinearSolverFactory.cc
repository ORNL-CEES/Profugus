//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testLinearSolverFactory.cc
 * \author Steven Hamilton
 * \brief  Test of LinearSolverFactory class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <time.h>

#include "comm/global.hh"
#include "../LinearSystem.hh"
#include "../LinearSystemFactory.hh"
#include "../LinearSolverFactory.hh"
#include "../RichardsonIteration.hh"
#include "../SyntheticAcceleration.hh"
//#include "../MonteCarloSolver.hh"
#include "../AdditiveSchwarzWrapper.hh"
#include "../BelosSolver.hh"
#include "../ChebyshevIteration.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace alea;

TEST(LinearSolverFactory, Basic)
{
    // Read ParameterList from file
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");

    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",10);

    pl->set("max_iterations",1000);
    pl->set("tolerance",1.0e-6);

    Teuchos::RCP<alea::LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();

    // Test Richardson
    Teuchos::RCP<OP> solver =
        alea::LinearSolverFactory::buildSolver("richardson",A,pl);

    Teuchos::RCP<RichardsonIteration> rich_it =
        Teuchos::rcp_dynamic_cast<RichardsonIteration>(solver);

    EXPECT_FALSE( rich_it == Teuchos::null );
    rich_it = Teuchos::null;

    /*
    // Test Monte Carlo
    solver = Teuchos::null;
    solver = alea::LinearSolverFactory::buildSolver("monte_carlo",A,pl);

    // Requesting monte_carlo actually returns an AdditiveSchwarzWrapper
    Teuchos::RCP<AdditiveSchwarzWrapper> mc_solver =
        Teuchos::rcp_dynamic_cast<AdditiveSchwarzWrapper>(solver);

    EXPECT_FALSE( mc_solver == Teuchos::null );
    mc_solver = Teuchos::null;

    // Test MCSA
    solver = Teuchos::null;
    solver = alea::LinearSolverFactory::buildSolver("mcsa",A,pl);

    Teuchos::RCP<SyntheticAcceleration> mcsa_solver =
        Teuchos::rcp_dynamic_cast<SyntheticAcceleration>(solver);

    EXPECT_FALSE( mcsa_solver == Teuchos::null );
    mcsa_solver = Teuchos::null;
    */

    // Test SyntheticAcceleration
    solver = Teuchos::null;
    pl->set("preconditioner","polynomial");
    solver = alea::LinearSolverFactory::buildSolver(
        "synthetic_acceleration",A,pl);

    Teuchos::RCP<SyntheticAcceleration> accel_solver =
        Teuchos::rcp_dynamic_cast<SyntheticAcceleration>(solver);

    EXPECT_FALSE( accel_solver == Teuchos::null );
    accel_solver = Teuchos::null;

    // Test Belos
    solver = Teuchos::null;
    solver = alea::LinearSolverFactory::buildSolver("belos",A,pl);

    Teuchos::RCP<BelosSolver> belos_solver =
        Teuchos::rcp_dynamic_cast<BelosSolver>(solver);

    EXPECT_FALSE( belos_solver == Teuchos::null );
    belos_solver = Teuchos::null;

    // Test Chebyshev
    solver = Teuchos::null;
    solver = alea::LinearSolverFactory::buildSolver("chebyshev",A,pl);

    Teuchos::RCP<ChebyshevIteration> cheby_solver =
        Teuchos::rcp_dynamic_cast<ChebyshevIteration>(solver);

    EXPECT_FALSE( cheby_solver == Teuchos::null );
}

