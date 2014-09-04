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
#include <Epetra_Vector.h>
#include <Thyra_EpetraLinearOp.hpp>
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

    void one_pe(const std::string &xmlfile)
    {
        if (nodes != 1)
            return;

        // database
        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        Teuchos::RCP<OP> A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);

        // wrap rhs into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_rhs = Teuchos::rcp(
            new Epetra_Vector(View, A->OperatorDomainMap(),&u1[0]) );

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),x) );
        int ret;
        ret = solver_op.Apply(*ep_rhs, *ep_x);
        EXPECT_EQ(0, ret);
        double tol = 1e-6;
        for( int i=0; i<4; ++i )
            EXPECT_SOFTEQ( x[i], sol[i], tol );

        // solve again and limit iterations
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        ret = solver_op.Apply(*ep_rhs,*ep_x);
        EXPECT_EQ(0, ret);
    }

    void four_pe(const std::string &xmlfile)
    {
        if (nodes != 4)
            return;

        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        Teuchos::RCP<OP> A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);

        // solve
        double x[1] = {0.0};
        double b[1] = {u1[node]};

        // wrap arrays into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_b = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),b) );
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),x) );

        int ret;
        ret = solver_op.Apply(*ep_b,*ep_x);
        EXPECT_EQ(0, ret);

        // reduce and check
        double global[4] = {0.0};
        global[node] = x[0];
        profugus::global_sum(global, 4);

        // check solution
        double tol = 1e-6;
        for( int i=0; i<4; ++i )
            EXPECT_SOFTEQ( global[i], sol[i], tol );
    }

    void one_pe_gen(const std::string &xmlfile)
    {
        if (nodes != 1)
            return;

        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        Teuchos::RCP<OP> A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);
        Teuchos::RCP<OP> B = linalg_traits::build_matrix<Matrix>("4x4_rhs",4);

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);
        solver_op.set_rhs_operator(B);

        // wrap rhs into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_rhs = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),&u2[0]) );

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),x) );
        int ret;
        ret = solver_op.Apply(*ep_rhs, *ep_x);
        EXPECT_EQ(0, ret);
        double tol = 1e-6;
        for( int i=0; i<4; ++i )
            EXPECT_SOFTEQ( x[i], sol[i], tol );

        // solve again and limit iterations
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        ret = solver_op.Apply(*ep_rhs,*ep_x);
        EXPECT_EQ(0, ret);
    }

    void four_pe_gen(const std::string &xmlfile)
    {
        if (nodes != 4)
            return;

        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        Teuchos::RCP<OP> A = linalg_traits::build_matrix<Matrix>("4x4_lhs",4);
        Teuchos::RCP<OP> B = linalg_traits::build_matrix<Matrix>("4x4_rhs",4);

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);
        solver_op.set_rhs_operator(B);

        // solve
        double x[1] = {0.0};
        double b[1] = {u2[node]};

        // wrap arrays into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_b = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),b) );
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),x) );

        int ret;
        ret = solver_op.Apply(*ep_b,*ep_x);
        EXPECT_EQ(0, ret);

        // reduce and check
        double global[4] = {0.0};
        global[node] = x[0];
        profugus::global_sum(global, 4);

        // check solution
        double tol = 1e-6;
        for( int i=0; i<4; ++i )
            EXPECT_SOFTEQ( global[i], sol[i], tol );
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
    one_pe("aztecoo.xml");
    four_pe("aztecoo.xml");
}

//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Belos)
{
    one_pe("belos.xml");
    four_pe("belos.xml");
}

//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Gen_Aztec)
{
    one_pe_gen("aztecoo.xml");
    four_pe_gen("aztecoo.xml");
}

//---------------------------------------------------------------------------//

TEST_F(Inverse_Operator_Test, Gen_Belos)
{
    one_pe_gen("belos.xml");
    four_pe_gen("belos.xml");
}

//---------------------------------------------------------------------------//
//                        end of tstInverseOperator.cc
//---------------------------------------------------------------------------//
