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

using namespace std;

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//
// Define a 4x4 matrix.

double A[][4] = {{10.0, 1.1, 2.0, 4.0},
                 { 1.1, 9.9, 2.1, 3.2},
                 { 0.8, 0.4, 5.3, 1.9},
                 { 0.3, 0.1, 0.4, 3.1}};

double B[][4] = {{0.56, 0.26, 0.51, 0.26},
                 {0.52, 0.13, 0.11, 0.41},
                 {0.73, 0.45, 0.40, 0.98},
                 {0.30, 0.44, 0.93, 0.35}};

double u1[] = {0.1, 0.3, 0.4, 0.9};

double u2[] = {  -1.233282375522885,
                -12.025836540161043,
                  5.290319979906284,
                  4.689586311471083};

double sol[] = {-0.102855551350840,
                -0.053521967514522,
                -0.013870314679620,
                 0.303792576783404};

//---------------------------------------------------------------------------//

class OperatorA : public Epetra_Operator
{
  public:
    typedef profugus::Decomposition::Comm Comm_t;
    typedef profugus::Decomposition::Map  Map_t;

  private:
    profugus::Decomposition d_map;
    int nodes, node;

  public:
    OperatorA(int num_elements)
        : d_map(num_elements)
        , nodes(profugus::nodes())
        , node(profugus::node())
    {

    }

    // Thyra interface
    int Apply(const Epetra_MultiVector &v, Epetra_MultiVector &y) const
    {
        if (nodes == 1)
        {
            for (int i = 0; i < 4; i++)
            {
                y[0][i] = 0.0;
                for (int j = 0; j < 4; j++)
                    y[0][i] += A[i][j] * v[0][j];
            }
        }

        if (nodes == 4)
        {
            // do a poor-man's gather
            double vv[4] = {0.0};
            vv[node] = v[0][0];
            profugus::global_sum(vv, 4);

            y[0][0] = 0.0;
            for (int j = 0; j < 4; j++)
                y[0][0] += A[node][j] * vv[j];
        }

        return 0;
    }

    const char* Label() const { return "A"; }
    const Comm_t& Comm() const { return d_map.comm(); }
    const Map_t& OperatorDomainMap() const { return d_map.map(); }
    const Map_t& OperatorRangeMap() const { return d_map.map(); }

    int SetUseTranspose(bool UseTranspose) { return 0; }
    int ApplyInverse(const Epetra_MultiVector &x,
                     Epetra_MultiVector &y) const { return 0; }
    double NormInf() const { return 0.0; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
};

//---------------------------------------------------------------------------//

class OperatorB : public Epetra_Operator
{
  public:
    typedef profugus::Decomposition::Comm Comm_t;
    typedef profugus::Decomposition::Map  Map_t;

  private:
    profugus::Decomposition d_map;
    int nodes, node;

  public:
    OperatorB(int num_elements)
        : d_map(num_elements)
        , nodes(profugus::nodes())
        , node(profugus::node())
    {

    }

    // Thyra interface
    int Apply(const Epetra_MultiVector &v, Epetra_MultiVector &y) const
    {
        if (nodes == 1)
        {
            for (int i = 0; i < 4; i++)
            {
                y[0][i] = 0.0;
                for (int j = 0; j < 4; j++)
                    y[0][i] += B[i][j] * v[0][j];
            }
        }

        if (nodes == 4)
        {
            // do a poor-man's gather
            double vv[4] = {0.0};
            vv[node] = v[0][0];
            profugus::global_sum(vv, 4);

            y[0][0] = 0.0;
            for (int j = 0; j < 4; j++)
                y[0][0] += B[node][j] * vv[j];
        }

        return 0;
    }

    const char* Label() const { return "A"; }
    const Comm_t& Comm() const { return d_map.comm(); }
    const Map_t& OperatorDomainMap() const { return d_map.map(); }
    const Map_t& OperatorRangeMap() const { return d_map.map(); }

    int SetUseTranspose(bool UseTranspose) { return 0; }
    int ApplyInverse(const Epetra_MultiVector &x,
                     Epetra_MultiVector &y) const { return 0; }
    double NormInf() const { return 0.0; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
};

//---------------------------------------------------------------------------//
// FIXTURES
//---------------------------------------------------------------------------//

class Inverse_Operator_Test : public testing::Test
{
  protected:
    typedef profugus::InverseOperator          InverseOperator;
    typedef InverseOperator::RCP_ParameterList RCP_ParameterList;

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
        Teuchos::RCP<OperatorA> A = Teuchos::rcp(new OperatorA(4));

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);

        // wrap rhs into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_rhs = Teuchos::rcp(
            new Epetra_Vector(View, A->OperatorDomainMap(),u1) );

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),x) );
        int ret;
        ret = solver_op.Apply(*ep_rhs, *ep_x);
        EXPECT_EQ(0, ret);
        EXPECT_TRUE(soft_equiv(x, x + 4, sol, sol + 4));

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
        Teuchos::RCP<OperatorA> A = Teuchos::rcp(new OperatorA(1));

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
        EXPECT_TRUE(soft_equiv(global, global + 4, sol, sol + 4));
    }

    void one_pe_gen(const std::string &xmlfile)
    {
        if (nodes != 1)
            return;

        db->set("linear_solver_xml_file", xmlfile);
        db->set("solver_type", std::string("Stratimikos"));

        // make the operator
        Teuchos::RCP<OperatorA> A = Teuchos::rcp(new OperatorA(4));
        Teuchos::RCP<OperatorB> B = Teuchos::rcp(new OperatorB(4));

        // make the solver
        InverseOperator solver_op(db);

        solver_op.set_operator(A);
        solver_op.set_rhs_operator(B);

        // wrap rhs into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_rhs = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),u2) );

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View,A->OperatorDomainMap(),x) );
        int ret;
        ret = solver_op.Apply(*ep_rhs, *ep_x);
        EXPECT_EQ(0, ret);
        EXPECT_TRUE(soft_equiv(x, x + 4, sol, sol + 4));

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
        Teuchos::RCP<OperatorA> A = Teuchos::rcp(new OperatorA(1));
        Teuchos::RCP<OperatorB> B = Teuchos::rcp(new OperatorB(1));

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
        EXPECT_TRUE(soft_equiv(global, global + 4, sol, sol + 4));
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
