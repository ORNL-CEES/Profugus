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
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

std::vector<double> rhs = {0.1, 0.3, 0.4, 0.9};

std::vector<double> sol = {-0.102855551350840, -0.053521967514522,
                           -0.013870314679620, 0.303792576783404};
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
    typedef Anasazi::MultiVecTraits<double,MV> MVT;
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

    void strat_test(const std::string &which,
                    const std::string &xmlfile)
    {
        // database
        d_db->set("linear_solver_xml_file", xmlfile);

        // make the operator
        d_A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);

        // make the solver
        Solver_t solver(d_db);

        solver.set_operator(d_A);

        // wrap rhs into Epetra MV
        Teuchos::RCP<MV> ep_rhs = linalg_traits::build_vector<MV>(4);
        linalg_traits::fill_vector<MV>(ep_rhs,rhs);

        std::vector<double> rhs_norm(1);
        MVT::MvNorm(*ep_rhs,rhs_norm);

        // solve
        Teuchos::RCP<MV> ep_x = linalg_traits::build_vector<MV>(4);
        std::vector<double> zero(4,0.0);
        linalg_traits::fill_vector<MV>(ep_x,zero);
        solver.set_tolerance(1.0e-8);
        solver.solve(ep_x, ep_rhs);
        linalg_traits::test_vector<MV>(ep_x,sol);
        EXPECT_TRUE( 10 >  solver.num_iters() );

        // solve again and limit iterations
        linalg_traits::fill_vector<MV>(ep_x,zero);
        solver.set_max_iters(1);
        solver.set_tolerance(1.0e-12);
        solver.solve(ep_x, ep_rhs);
        EXPECT_TRUE( 10 >  solver.num_iters() );

        // solve again and limit tolerance
        linalg_traits::fill_vector<MV>(ep_x,zero);
        solver.set_max_iters(1000);
        solver.set_tolerance(0.1);
        solver.solve(ep_x, ep_rhs);
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
    strat_test("AztecOO", "aztecoo.xml");
}

//---------------------------------------------------------------------------//

TEST_F(StratimikosSolver_Test, Belos)
{
    strat_test("Belos", "belos.xml");
}

//---------------------------------------------------------------------------//
#ifdef USE_MCLS
TEST_F(StratimikosSolver_Test, MCLS)
{
    strat_test("MCLS", "mcls.xml");
}
#endif

//---------------------------------------------------------------------------//
//                 end of tstStratimikosSolver.cc
//---------------------------------------------------------------------------//
