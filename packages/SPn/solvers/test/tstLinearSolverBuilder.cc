//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstLinearSolverBuilder.cc
 * \author Steven Hamilton
 * \brief  LinearSolverBuilder unit test.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <SPn/config.h>

#include "../LinearSolverBuilder.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class SolverBuilderTest : public testing::Test
{
  protected:

    typedef profugus::LinearSolverBuilder Builder;
    typedef Builder::RCP_ParameterList    RCP_ParameterList;
    typedef Builder::RCP_LinearSolver     RCP_LinearSolver;
    typedef Builder::MV                   MV;
    typedef Builder::OP                   OP;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
    }

    void build_solver( RCP_ParameterList db )
    {
        d_solver = Builder::build_solver(db);
    }

  protected:
    int node;
    int nodes;

    RCP_LinearSolver d_solver;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST_F(SolverBuilderTest, basic)
{
    RCP_ParameterList db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));

    // Default solver is Richardson
    build_solver(db);
    EXPECT_EQ("Exnihilo Richardson", d_solver->solver_label());

    //
    // Exnihilo solver by specifying exnihilo_solver
    //

    db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));
    db->set("exnihilo_solver", std::string("Richardson"));
    build_solver(db);
    EXPECT_EQ("Exnihilo Richardson", d_solver->solver_label());

    //
    // Exnihilo solver by specifying solver_type
    //

    db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));
    db->set("solver_type", std::string("Exnihilo"));
    build_solver(db);
    EXPECT_EQ("Exnihilo Richardson", d_solver->solver_label());

    //
    // Stratimikos solver (default is AztecOO)
    //

    db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));
    db->set("solver_type", std::string("Stratimikos"));
    build_solver(db);
    EXPECT_EQ("Stratimikos AztecOO", d_solver->solver_label());
}

//---------------------------------------------------------------------------//
//                        end of tstLinearSolverBuilder.cc
//---------------------------------------------------------------------------//
