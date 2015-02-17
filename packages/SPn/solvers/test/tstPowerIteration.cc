//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstPowerIteration.cc
 * \author Steven Hamilton
 * \brief  PowerIteration unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <SPn/config.h>

#include "../PowerIteration.hh"
#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

template <class T>
class PowerIterationTest : public ::testing::Test
{
  protected:

    typedef typename T::MV              MV;
    typedef typename T::MATRIX          MATRIX;
    typedef profugus::PowerIteration<T> PowerIteration;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Parallelism
        node  = profugus::node();
        nodes = profugus::nodes();

        // Build matrix
        d_N = 20;
        d_A = linalg_traits::build_matrix<T>("laplacian",d_N);

        // Build eigenvector
        d_x = linalg_traits::build_vector<T>(d_N);
        std::vector<double> one(d_N,1.0);
        linalg_traits::fill_vector<T>(d_x,one);

        // Create options database
        d_db = Teuchos::rcp(new Teuchos::ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);

        // Build solver
        d_solver = Teuchos::rcp(new PowerIteration(d_db));
        CHECK(!d_solver.is_null());
        d_solver->set_operator(d_A);
    }

    void solve()
    {
        d_solver->solve(d_lambda,d_x);
        d_iters = d_solver->num_iters();
        d_converged = d_solver->converged();
    }

  protected:
    int node;
    int nodes;
    int d_N;

    Teuchos::RCP<Teuchos::ParameterList> d_db;
    Teuchos::RCP<MATRIX>         d_A;
    Teuchos::RCP<MV>             d_x;
    Teuchos::RCP<PowerIteration> d_solver;
    double                       d_lambda;

    int d_iters;
    bool d_converged;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(PowerIterationTest, MyTypes);

TYPED_TEST(PowerIterationTest, basic)
{
    // Run two iterations and stop
    this->solve();

    // Make sure solver reports that two iterations were performed
    //  and that it is not converged
    EXPECT_EQ( 2, this->d_iters );
    EXPECT_TRUE( !this->d_converged );

    // Reset initial vector and re-solve
    this->d_solver->set_max_iters(1000);
    std::vector<double> one(this->d_N,1.0);
    linalg_traits::fill_vector<TypeParam>(this->d_x,one);
    this->solve();

    EXPECT_EQ( 261, this->d_iters );
    EXPECT_TRUE( this->d_converged );
    EXPECT_SOFTEQ( this->d_lambda, 3.911145611572282, 1.0e-6 );

    // Compare against reference solution from Matlab
    std::vector<double> ref =
            {9.096342209328087e-02,  -1.738443448352310e-01,
             2.412784344502512e-01,  -2.872738761376646e-01,
             3.077437731144456e-01,  -3.008692853090870e-01,
             2.672612412471236e-01,  -2.099058632158057e-01,
             1.338994282792137e-01,  -4.599544168741217e-02,
             -4.599544168741519e-02,  1.338994282792167e-01,
             -2.099058632158085e-01,  2.672612412471261e-01,
             -3.008692853090894e-01,  3.077437731144476e-01,
             -2.872738761376664e-01,  2.412784344502525e-01,
             -1.738443448352319e-01,  9.096342209328130e-02};

    linalg_traits::test_vector<TypeParam>(this->d_x,ref);

    // Solve again, should return without iterating
    this->solve();

    EXPECT_EQ( 1, this->d_iters );
    EXPECT_TRUE( this->d_converged );

    // Make sure solution didn't change
    linalg_traits::test_vector<TypeParam>(this->d_x,ref);
}

//---------------------------------------------------------------------------//
//                        end of tstPowerIteration.cc
//---------------------------------------------------------------------------//
