//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstEigenvalueSolverBuilder.cc
 * \author Steven Hamilton
 * \brief  EigenvalueSolverBuilder unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include <SPn/config.h>

#include "../EigenvalueSolverBuilder.hh"
#include "../Arnoldi.hh"
#include "../Davidson_Eigensolver.hh"
#include "../PowerIteration.hh"
#include "../RayleighQuotient.hh"
#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class EigenvalueSolverBuilderTest : public ::testing::Test
{
  protected:

    typedef Epetra_MultiVector                       MV;
    typedef Epetra_Operator                          OP;
    typedef profugus::EigenvalueSolverBuilder<MV,OP> Builder;
    typedef Builder::RCP_ParameterList               RCP_ParameterList;
    typedef Builder::RCP_EigenvalueSolver            RCP_EigenvalueSolver;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        int num_global = 4;
        d_A = linalg_traits::build_matrix<Epetra_CrsMatrix>("laplacian",num_global);
        d_B = linalg_traits::build_matrix<Epetra_CrsMatrix>("diagonal",num_global);
    }

  protected:

    RCP_EigenvalueSolver           d_solver;
    Teuchos::RCP<Epetra_CrsMatrix> d_A;
    Teuchos::RCP<Epetra_CrsMatrix> d_B;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST_F(EigenvalueSolverBuilderTest, basic)
{
    RCP_ParameterList db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));

    // Default standard eigenvalue solver is Arnoldi
    d_solver = Builder::build_solver(db,d_A);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());
    Teuchos::RCP<profugus::Arnoldi<MV,OP> > arnoldi =
        Teuchos::rcp_dynamic_cast<profugus::Arnoldi<MV,OP> >(d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    // Default generalized eigenvalue solver is Arnoldi (for now)
    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());
    arnoldi = Teuchos::rcp_dynamic_cast<profugus::Arnoldi<MV,OP> >(d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    // Make sure "Arnoldi" keyword is recognized by both functions
    db->set("eigensolver",std::string("Arnoldi"));

    d_solver = Builder::build_solver(db,d_A);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());
    arnoldi = Teuchos::rcp_dynamic_cast<profugus::Arnoldi<MV,OP> >(d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());
    arnoldi = Teuchos::rcp_dynamic_cast<profugus::Arnoldi<MV,OP> >(d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    // Power iteration for both standard and generalized problems
    db->set("eigensolver",std::string("Power"));

    d_solver = Builder::build_solver(db,d_A);
    EXPECT_EQ("Power Iteration",d_solver->solver_label());
    Teuchos::RCP<profugus::PowerIteration<MV,OP> > power =
        Teuchos::rcp_dynamic_cast<profugus::PowerIteration<MV,OP> >(d_solver);
    EXPECT_TRUE( power != Teuchos::null );

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Power Iteration",d_solver->solver_label());
    power = Teuchos::rcp_dynamic_cast<profugus::PowerIteration<MV,OP> >(d_solver);
    EXPECT_TRUE( power != Teuchos::null );

    // Rayleigh quotient iteration for generalized problem
    db->set("eigensolver",std::string("RQI"));

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Rayleigh Quotient",d_solver->solver_label());
    Teuchos::RCP<profugus::RayleighQuotient<MV,OP> > rqi =
        Teuchos::rcp_dynamic_cast<profugus::RayleighQuotient<MV,OP> >(d_solver);
    EXPECT_TRUE( rqi != Teuchos::null );

    // Davidson for generalized problem
    db->set("eigensolver",std::string("Davidson"));

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Davidson",d_solver->solver_label());
    Teuchos::RCP<profugus::Davidson_Eigensolver<MV,OP> > davidson =
        Teuchos::rcp_dynamic_cast<profugus::Davidson_Eigensolver<MV,OP> >(d_solver);
    EXPECT_TRUE( davidson!= Teuchos::null );

}

//---------------------------------------------------------------------------//
//                        end of tstEigenvalueSolverBuilder.cc
//---------------------------------------------------------------------------//
