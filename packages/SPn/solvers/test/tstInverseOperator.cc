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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include "../Decomposition.hh"
#include "../InverseOperator.hh"
#include "LinAlgTraits.hh"

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

template <class T>
class Inverse_Operator_Test : public testing::Test
{
  protected:

    typedef typename linalg_traits::traits_types<T>::MV     MV;
    typedef typename linalg_traits::traits_types<T>::OP     OP;
    typedef typename linalg_traits::traits_types<T>::Matrix Matrix;
    typedef profugus::InverseOperator<MV,OP>                InverseOperator;
    typedef Anasazi::OperatorTraits<double,MV,OP>           OPT;

  protected:

    void SetUp()
    {
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
        OPT::Apply(solver_op,*ep_rhs,*ep_x);
        linalg_traits::test_vector<MV>(ep_x,sol);

        // solve again and limit iterations
        std::vector<double> zero(N,0.0);
        linalg_traits::fill_vector<MV>(ep_x,zero);
        OPT::Apply(solver_op,*ep_rhs,*ep_x);
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
        OPT::Apply(solver_op,*ep_rhs,*ep_x);
        linalg_traits::test_vector(ep_x,sol);

        // solve again and limit iterations
        std::vector<double> zero(N,0.0);
        linalg_traits::fill_vector<MV>(ep_x,zero);
        OPT::Apply(solver_op,*ep_rhs,*ep_x);
    }

  protected:

    Teuchos::RCP<Teuchos::ParameterList> db;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
typedef ::testing::Types<Epetra_MultiVector,Tpetra_MultiVector> MyTypes;
TYPED_TEST_CASE(Inverse_Operator_Test, MyTypes);

TYPED_TEST(Inverse_Operator_Test, Aztec)
{
    this->std_test("aztecoo.xml");
}

//---------------------------------------------------------------------------//

TYPED_TEST(Inverse_Operator_Test, Belos)
{
    this->std_test("belos.xml");
}

//---------------------------------------------------------------------------//

TYPED_TEST(Inverse_Operator_Test, Gen_Aztec)
{
    this->gen_test("aztecoo.xml");
}

//---------------------------------------------------------------------------//

TYPED_TEST(Inverse_Operator_Test, Gen_Belos)
{
    this->gen_test("belos.xml");
}

//---------------------------------------------------------------------------//
//                        end of tstInverseOperator.cc
//---------------------------------------------------------------------------//
