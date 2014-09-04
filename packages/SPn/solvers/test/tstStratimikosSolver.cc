//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstStratimikosSolver.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 17 21:12:05 2014
 * \brief  StratimikosSolver Unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <string>

#include "../Decomposition.hh"
#include "../StratimikosSolver.hh"

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

double rhs[] = {0.1, 0.3, 0.4, 0.9};

double sol[] = {-0.102855551350840,
                -0.053521967514522,
                -0.013870314679620,
                0.303792576783404};

int nodes, node;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class StratimikosSolver_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef Epetra_MultiVector                 MV;
    typedef Epetra_Operator                    OP;
    typedef Epetra_CrsMatrix                   Matrix;
    typedef profugus::StratimikosSolver<MV,OP> Solver_t;
    typedef Solver_t::ParameterList            ParameterList;
    typedef Solver_t::RCP_ParameterList        RCP_ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        d_db = Teuchos::rcp(new ParameterList("test"));
    }

    void one_pe(const std::string &which,
                const std::string &xmlfile)
    {
        if (nodes != 1)
        {
            SUCCEED() << "Test set for 1 processor only.";
            return;
        }

        // database
        d_db->set("linear_solver_xml_file", xmlfile);

        // make the operator
        d_A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);

        // make the solver
        Solver_t solver(d_db);

        solver.set_operator(d_A);

        // wrap rhs into Epetra MV
        Teuchos::RCP<MV> ep_rhs = linalg_traits::build_vector<MV>(4);
        for( int i=0; i<4; ++i )
            ep_rhs->ReplaceGlobalValue(i,0,rhs[i]);

        std::vector<double> rhs_norm(1);
        ep_rhs->Norm2(&rhs_norm[0]);

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<MV> ep_x = linalg_traits::build_vector<MV>(4);
        for( int i=0; i<4; ++i )
            ep_x->ReplaceGlobalValue(i,0,x[i]);
        solver.set_tolerance(1.0e-8);
        solver.solve(ep_x, ep_rhs);
        EXPECT_TRUE(soft_equiv((*ep_x)[0], (*ep_x)[0] + 4, sol, sol + 4, 1.0e-6));
        EXPECT_TRUE( 10 >  solver.num_iters() );

        // solve again and limit iterations
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        solver.set_max_iters(1);
        solver.set_tolerance(1.0e-12);
        solver.solve(ep_x, ep_rhs);
        EXPECT_TRUE( 10 >  solver.num_iters() );

        // solve again and limit tolerance
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        solver.set_max_iters(1000);
        solver.set_tolerance(0.1);
        solver.solve(ep_x, ep_rhs);
        EXPECT_TRUE( 10 >  solver.num_iters() );
    }

    void four_pe(const std::string &which,
                 const std::string &xmlfile)
    {
        if (nodes != 4)
        {
            SUCCEED() << "Test set for 4 processors only.";
            return;
        }

        // database
        d_db->set("linear_solver_xml_file", xmlfile);

        // make the operator
        d_A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);

        // make the solver
        Solver_t solver(d_db);

        solver.set_operator(d_A);

        // wrap arrays into Epetra MV
        Teuchos::RCP<MV> ep_b = linalg_traits::build_vector<MV>(4);
        for( int i=0; i<2; ++i )
        {
            int global = d_A->GRID(i);
            ep_b->ReplaceGlobalValue(global,0,rhs[global]);
        }
        Teuchos::RCP<MV> ep_x = linalg_traits::build_vector<MV>(4);
        ep_x->PutScalar(0.0);

        solver.solve(ep_x, ep_b);

        // reduce and check
        double global[4] = {0.0};
        global[node] = (*ep_x)[0][0];
        profugus::global_sum(global, 4);

        // check solution
        EXPECT_TRUE(soft_equiv(global, global + 4, sol, sol + 4, 1.0e-5));
        EXPECT_TRUE( 10 >  solver.num_iters() );
    }

  protected:
    // >>> Data that get re-initialized between tests

    RCP_ParameterList d_db;
    Teuchos::RCP<Matrix> d_A;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(StratimikosSolver_Test, Aztec)
{
    one_pe("AztecOO", "aztecoo.xml");
    four_pe("AztecOO", "aztecoo.xml");
}

//---------------------------------------------------------------------------//

TEST_F(StratimikosSolver_Test, Belos)
{
    one_pe("Belos", "belos.xml");
    four_pe("Belos", "belos.xml");
}

//---------------------------------------------------------------------------//
#ifdef USE_MCLS
TEST_F(StratimikosSolver_Test, MCLS)
{
    one_pe("MCLS", "mcls.xml");
    four_pe("MCLS", "mcls.xml");
}
#endif

//---------------------------------------------------------------------------//
//                 end of tstStratimikosSolver.cc
//---------------------------------------------------------------------------//
