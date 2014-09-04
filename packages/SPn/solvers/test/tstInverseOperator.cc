//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstInverseOperator.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Jan 27 11:56:51 2009
 * \brief  InverseOperator class test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include "../Decomposition.hh"
#include "../InverseOperator.hh"
#include "LinAlgTraits.hh"

using namespace std;

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

std::vector<double> u1 = {0.1, 0.3, 0.4, 0.9};

std::vector<double> u2 = {  -1.233282375522885,
                           -12.025836540161043,
                             5.290319979906284,
                             4.689586311471083};

std::vector<double> sol = {-0.102855551350840,
                           -0.053521967514522,
                           -0.013870314679620,
                            0.303792576783404};

//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
// FIXTURES
//---------------------------------------------------------------------------//

class Inverse_Operator_Test : public testing::Test
{
  protected:
    typedef Epetra_MultiVector                   MV;
    typedef Epetra_Operator                      OP;
    typedef Epetra_CrsMatrix                     Matrix;
    typedef profugus::InverseOperator<MV,OP>     InverseOperator;
    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_ParameterList;

  protected:

    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
        db    = Teuchos::rcp(new Teuchos::ParameterList("test"));
    }

    void std_test(const std::string &xmlfile)
    {
        // database
        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        int N = 4;
        Teuchos::RCP<OP> A = linalg_traits::build_matrix<Matrix>("4x4_lhs",N);

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);

        // wrap rhs into MV
        Teuchos::RCP<MV> ep_rhs = linalg_traits::build_vector<MV>(N);
        linalg_traits::fill_vector<MV>(ep_rhs,u1);

        // solve
        Teuchos::RCP<MV> ep_x = linalg_traits::build_vector<MV>(N);
        int ret;
        ret = solver_op.Apply(*ep_rhs, *ep_x);
        EXPECT_EQ(0, ret);
        linalg_traits::test_vector<MV>(ep_x,sol);

        // solve again and limit iterations
        std::vector<double> zero(N,0.0);
        linalg_traits::fill_vector<MV>(ep_x,zero);
        ret = solver_op.Apply(*ep_rhs,*ep_x);
        EXPECT_EQ(0, ret);
    }

    void gen_test(const std::string &xmlfile)
    {
        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        int N = 4;
        Teuchos::RCP<OP> A = linalg_traits::build_matrix<Matrix>("4x4_lhs",N);
        Teuchos::RCP<OP> B = linalg_traits::build_matrix<Matrix>("4x4_rhs",N);

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);
        solver_op.set_rhs_operator(B);

        // wrap rhs into MV
        Teuchos::RCP<MV> ep_rhs = linalg_traits::build_vector<MV>(N);
        linalg_traits::fill_vector<MV>(ep_rhs,u2);

        // solve
        Teuchos::RCP<MV> ep_x = linalg_traits::build_vector<MV>(N);
        int ret;
        ret = solver_op.Apply(*ep_rhs, *ep_x);
        EXPECT_EQ(0, ret);
        linalg_traits::test_vector(ep_x,sol);

        // solve again and limit iterations
        std::vector<double> zero(N,0.0);
        linalg_traits::fill_vector<MV>(ep_x,zero);
        ret = solver_op.Apply(*ep_rhs,*ep_x);
        EXPECT_EQ(0, ret);
    }

  protected:
    int nodes, node;

    RCP_ParameterList db;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Aztec)
{
    std_test("aztecoo.xml");
}

//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Belos)
{
    std_test("belos.xml");
}

//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Gen_Aztec)
{
    gen_test("aztecoo.xml");
}

//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Gen_Belos)
{
    gen_test("belos.xml");
}

//---------------------------------------------------------------------------//
//                        end of tstInverseOperator.cc
//---------------------------------------------------------------------------//
