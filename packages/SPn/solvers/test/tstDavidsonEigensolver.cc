//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstDavidsonEigensolver.cc
 * \author Steven Hamilton
 * \brief  Davidson_Eigensolver unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>
#include <string>

#include <SPn/config.h>

#include "../Davidson_Eigensolver.hh"
#include "LinAlgTraits.hh"

// Reference solution from matlab
double ref_eigenvalue = 0.4890748754542557;

std::vector<double> ref_eigenvector = {
    4.599544191e-02, 9.096342166e-02, 1.338994289e-01, 1.738443441e-01,
    2.099058640e-01, 2.412784337e-01, 2.672612419e-01, 2.872738756e-01,
    3.008692856e-01, 3.077437730e-01, 3.077437730e-01, 3.008692856e-01,
    2.872738756e-01, 2.672612419e-01, 2.412784337e-01, 2.099058640e-01,
    1.738443441e-01, 1.338994289e-01, 9.096342166e-02, 4.599544191e-02 };

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

template <class T>
class DavidsonTest : public ::testing::Test
{
  protected:

    typedef typename T::MV       MV;
    typedef typename T::OP       OP;
    typedef typename T::MATRIX   MATRIX;

    typedef profugus::Davidson_Eigensolver<T> Davidson_Eigensolver;

  protected:

    // Initialization that are performed for each test
    void build_solver()
    {
        using Teuchos::rcp;

        // Build an map
        d_N = 20;
        d_A = linalg_traits::build_matrix<T>("shifted_laplacian",d_N);
        d_B = linalg_traits::build_matrix<T>("scaled_identity",d_N);

        // Build eigenvector
        d_x = linalg_traits::build_vector<T>(d_N);

        // Create options database
        d_db = rcp(new ParameterList("test"));
        d_db->set("tolerance",1e-8);

        // Build solver
        d_solver = rcp(new Davidson_Eigensolver(d_db,d_A,d_B));
        CHECK(!d_solver.is_null());
    }

    void solve()
    {
        d_solver->solve(d_lambda,d_x);
        d_iters = d_solver->num_iters();
        d_converged = d_solver->converged();
    }

  protected:
    int d_N;

    Teuchos::RCP<Teuchos::ParameterList> d_db;
    Teuchos::RCP<MATRIX>                 d_A;
    Teuchos::RCP<MATRIX>                 d_B;
    Teuchos::RCP<MV>                     d_x;
    Teuchos::RCP<Davidson_Eigensolver>   d_solver;
    double                               d_lambda;

    int d_iters;
    bool d_converged;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(DavidsonTest, MyTypes);

TYPED_TEST(DavidsonTest, basic)
{
    double eig_tol = 1e-10;
    double vec_tol = 1e-7;
    this->build_solver();

    // Solve and check convergence
    std::vector<double> one(this->d_N,1.0);
    linalg_traits::fill_vector<TypeParam>(this->d_x,one);
    this->d_lambda = 1.0;
    this->solve();

    EXPECT_EQ( 9, this->d_iters ); // Heuristic
    EXPECT_TRUE( this->d_converged );
    EXPECT_SOFTEQ( this->d_lambda, ref_eigenvalue, eig_tol );

    linalg_traits::set_sign<TypeParam>(this->d_x);
    linalg_traits::test_vector<TypeParam>(this->d_x,ref_eigenvector);

    // Solve again, should return without iterating
    this->solve();

    EXPECT_EQ( 0, this->d_iters );
    EXPECT_TRUE( this->d_converged );
    EXPECT_SOFTEQ( this->d_lambda, ref_eigenvalue, eig_tol );

    // Make sure solution didn't change
    linalg_traits::set_sign<TypeParam>(this->d_x);
    linalg_traits::test_vector<TypeParam>(this->d_x,ref_eigenvector);
}

//---------------------------------------------------------------------------//
//                 end of tstDavidsonEigensolver.cc
//---------------------------------------------------------------------------//
