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
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

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

class Operator
{
  public:
    typedef profugus::Decomposition Decomposition;

  private:
    Decomposition d_map;
    Teuchos::RCP<Epetra_CrsMatrix> d_crs_matrix;

  public:
    Operator(int num_elements)
        : d_map(num_elements)
    {
	int num_cols = 4;
	d_crs_matrix = Teuchos::rcp(
	    new Epetra_CrsMatrix(Copy, d_map.map(), num_cols) );
	Teuchos::Array<int> indices( num_cols );
	for ( int i = 0; i < num_cols; ++i )
	{
	    indices[i] = i;
	}
	if ( 4 == num_elements )
	{
	    for ( int i = 0; i < num_elements; ++i )
	    {
		d_crs_matrix->InsertGlobalValues( i, num_cols, matrix[i], indices.getRawPtr() );
	    }
	}
	else if ( 1 == num_elements )
	{
	    int my_rank = d_map.comm().MyPID();
	    d_crs_matrix->InsertGlobalValues( 
		my_rank , num_cols, &matrix[my_rank][0], indices.getRawPtr() );
	}
	d_crs_matrix->FillComplete();
    }

    Teuchos::RCP<Epetra_CrsMatrix> getOperator() const
    { return d_crs_matrix; }
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

        solver.set_operator(A->getOperator());

        // wrap rhs into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_rhs = Teuchos::rcp(
            new Epetra_Vector(View, A->getOperator()->OperatorDomainMap(), rhs));

        // solve
        double x[4] = {0.0};
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View, A->getOperator()->OperatorDomainMap(), x));
        solver.set_tolerance(1.0e-8);
        solver.solve(ep_x, ep_rhs);
        EXPECT_TRUE(soft_equiv(x, x + 4, sol, sol + 4, 1.0e-6));
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
        db->set("linear_solver_xml_file", xmlfile);

        // make the operator
        A = Teuchos::rcp(new Operator(1));

        // make the solver
        Solver_t solver(db);

        solver.set_operator(A->getOperator());

        // solve
        double x[1] = {0.0};
        double b[1] = {rhs[node]};

        // wrap arrays into Epetra MV
        Teuchos::RCP<Epetra_Vector> ep_b = Teuchos::rcp(
            new Epetra_Vector(View, A->getOperator()->OperatorDomainMap(), b));
        Teuchos::RCP<Epetra_Vector> ep_x = Teuchos::rcp(
            new Epetra_Vector(View, A->getOperator()->OperatorDomainMap(), x));

        solver.solve(ep_x, ep_b);

        // reduce and check
        double global[4] = {0.0};
        global[node] = x[0];
        profugus::global_sum(global, 4);

        // check solution
        EXPECT_TRUE(soft_equiv(global, global + 4, sol, sol + 4, 1.0e-5));
        EXPECT_TRUE( 10 >  solver.num_iters() );
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
