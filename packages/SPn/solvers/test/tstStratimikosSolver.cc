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

#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//
// Define a 4x4 matrix.

double matrix[][4] = {{10.0, 1.1, 2.0, 4.0},
                      { 1.1, 9.9, 2.1, 3.2},
                      { 0.8, 0.4, 5.3, 1.9},
                      { 0.3, 0.1, 0.4, 3.1}};

double rhs[] = {0.1, 0.3, 0.4, 0.9};

double sol[] = {-0.102855551350840,
                -0.053521967514522,
                -0.013870314679620,
                0.303792576783404};

int nodes, node;

//---------------------------------------------------------------------------//

class Operator : public Epetra_Operator
{
  public:
    typedef profugus::Decomposition Decomposition;

  private:
    Decomposition d_map;

  public:
    Operator(int num_elements)
        : d_map(num_elements)
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
                    y[0][i] += matrix[i][j] * v[0][j];
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
                y[0][0] += matrix[node][j] * vv[j];
        }

        return 0;
    }

    const char* Label() const { return "matrix"; }
    const Decomposition::Comm& Comm() const { return d_map.comm(); }
    const Decomposition::Map& OperatorDomainMap() const { return d_map.map(); }
    const Decomposition::Map& OperatorRangeMap() const { return d_map.map(); }

    int SetUseTranspose(bool UseTranspose) { return 0; }
    int ApplyInverse(const Epetra_MultiVector &x,
                     Epetra_MultiVector &y) const { return 0; }
    double NormInf() const { return 0.0; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class StratimikosSolver_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::StratimikosSolver Solver_t;
    typedef Teuchos::RCP<Operator>      RCP_Operator;
    typedef Solver_t::ParameterList     ParameterList;
    typedef Solver_t::RCP_ParameterList RCP_ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        db = Teuchos::rcp(new ParameterList("test"));
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
        db->set("linear_solver_xml_file", xmlfile);

        // make the operator
        A = Teuchos::rcp(new Operator(4));

        // make the solver
        Solver_t solver(db);

        solver.set_operator(A);

        // wrap rhs into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_rhs = Teuchos::rcp(
            new Epetra_Vector(View, A->OperatorDomainMap(), rhs));

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View, A->OperatorDomainMap(), x));
        solver.solve(ep_x, ep_rhs);
        EXPECT_TRUE(soft_equiv(x, x + 4, sol, sol + 4));
        EXPECT_EQ(4, solver.num_iters());

        // solve again and limit iterations
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        solver.set_max_iters(1);
        solver.set_tolerance(1.0e-12);
        solver.solve(ep_x, ep_rhs);
        EXPECT_EQ(1, solver.num_iters());

        // solve again and limit tolerance
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        solver.set_max_iters(1000);
        solver.set_tolerance(0.1);
        solver.solve(ep_x, ep_rhs);
        EXPECT_EQ(3, solver.num_iters());
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
        db->set("linear_solver_xml_file", xmlfile);

        // make the operator
        A = Teuchos::rcp(new Operator(1));

        // make the solver
        Solver_t solver(db);

        solver.set_operator(A);

        // solve
        double x[1] = {0.0};
        double b[1] = {rhs[node]};

        // wrap arrays into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_b = Teuchos::rcp(
            new Epetra_Vector(View, A->OperatorDomainMap(), b));
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View, A->OperatorDomainMap(), x));

        solver.solve(ep_x, ep_b);

        // reduce and check
        double global[4] = {0.0};
        global[node] = x[0];
        profugus::global_sum(global, 4);

        // check solution
        EXPECT_TRUE(soft_equiv(global, global + 4, sol, sol + 4));
    }

  protected:
    // >>> Data that get re-initialized between tests

    RCP_ParameterList db;
    RCP_Operator      A;
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
//                 end of tstStratimikosSolver.cc
//---------------------------------------------------------------------------//
