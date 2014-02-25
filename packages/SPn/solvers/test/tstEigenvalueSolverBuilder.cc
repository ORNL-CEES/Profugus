//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstEigenvalueSolverBuilder.cc
 * \author Steven Hamilton
 * \brief  EigenvalueSolverBuilder unit-test.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
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

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class EigenvalueSolverBuilderTest : public ::testing::Test
{
  protected:

    typedef profugus::EigenvalueSolverBuilder Builder;
    typedef Builder::RCP_ParameterList        RCP_ParameterList;
    typedef Builder::RCP_EigenvalueSolver     RCP_EigenvalueSolver;
    typedef Builder::MV                       MV;
    typedef Builder::OP                       OP;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Build Epetra communicator
#ifdef COMM_MPI
        Epetra_MpiComm comm(profugus::communicator);
#else
        Epetra_SerialComm comm;
#endif
        int num_global = profugus::nodes();
        Epetra_Map map(num_global,0,comm);
        d_A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,map,1) );
        d_B = Teuchos::rcp( new Epetra_CrsMatrix(Copy,map,1) );
        std::vector<int> ind(1);
        std::vector<double> val(1);
        ind[0] = profugus::node();
        val[0] = static_cast<double>((profugus::node()+1)*2);
        d_A->InsertMyValues(0,1,&val[0],&ind[0]);
        d_A->FillComplete();
        val[0] = 2.0;
        d_B->InsertMyValues(0,1,&val[0],&ind[0]);
        d_B->FillComplete();
    }

  protected:
    int node;
    int nodes;

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

    // Default generalized eigenvalue solver is Arnoldi (for now)
    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());

    // Make sure "Arnoldi" keyword is recognized by both functions
    db->set("eigensolver",std::string("Arnoldi"));

    d_solver = Builder::build_solver(db,d_A);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Arnoldi",d_solver->solver_label());

    // Power iteration for both standard and generalized problems
    db->set("eigensolver",std::string("Power"));

    d_solver = Builder::build_solver(db,d_A);
    EXPECT_EQ("Power Iteration",d_solver->solver_label());

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Power Iteration",d_solver->solver_label());

    // Rayleigh quotient iteration for generalized problem
    db->set("eigensolver",std::string("RQI"));

    d_solver = Builder::build_solver(db,d_A,d_B);
    EXPECT_EQ("Rayleigh Quotient",d_solver->solver_label());
}

//---------------------------------------------------------------------------//
//                        end of tstEigenvalueSolverBuilder.cc
//---------------------------------------------------------------------------//
