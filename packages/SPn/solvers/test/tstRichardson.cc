//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstRichardson.cc
 * \author Steven Hamilton
 * \brief  Richardson unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include <SPn/config.h>
#include "../Richardson.hh"

#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

template <class T>
class RichardsonTest : public testing::Test
{
  protected:

    typedef typename linalg_traits::traits_types<T>::MV       MV;
    typedef typename linalg_traits::traits_types<T>::OP       OP;
    typedef typename linalg_traits::traits_types<T>::Matrix   Matrix;

    typedef profugus::Richardson<MV,OP>   Richardson;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Build a map
        d_N = 8;
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",d_N);

        // Build lhs and rhs vectors
        d_x = linalg_traits::build_vector<MV>(d_N);
        d_b = linalg_traits::build_vector<MV>(d_N);
        std::vector<double> vals(d_N);
        for( int i=0; i<d_N; ++i )
            vals[i] = static_cast<double>(4*(8-i));
        linalg_traits::fill_vector<MV>(d_b,vals);

        // Create options database
        d_db = Teuchos::rcp(new Teuchos::ParameterList("test"));
        d_db->set("tolerance",1e-8);
        d_db->set("max_itr",2);
        d_db->set("Damping Factor",0.5);

        // Build solver
        d_solver = Teuchos::rcp(new Richardson(d_db));
        CHECK(!d_solver.is_null());
        d_solver->set_operator(d_A);
    }

    void solve()
    {
        d_solver->solve(d_x,d_b);
        d_iters = d_solver->num_iters();
        d_converged = d_solver->converged();
    }

  protected:
    int d_N;

    Teuchos::RCP<Teuchos::ParameterList> d_db;
    Teuchos::RCP<Matrix>                 d_A;
    Teuchos::RCP<MV>                     d_x;
    Teuchos::RCP<MV>                     d_b;
    Teuchos::RCP<Richardson>             d_solver;

    int d_iters;
    bool d_converged;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

typedef ::testing::Types<Epetra_MultiVector,Tpetra_MultiVector> MyTypes;
TYPED_TEST_CASE(RichardsonTest, MyTypes);

TYPED_TEST(RichardsonTest, basic)
{
    typedef typename TestFixture::MV       MV;

    // Run two iterations and stop
    this->solve();

    // Make sure solver reports that two iterations were performed
    //  and that it is not converged
    EXPECT_EQ( 2, this->d_iters );
    EXPECT_TRUE( !this->d_converged );

    // Reset initial vector and re-solve
    this->d_solver->set_max_iters(1000);
    std::vector<double> one(this->d_N,1.0);
    linalg_traits::fill_vector<MV>(this->d_x,one);
    this->solve();

    EXPECT_EQ( 294, this->d_iters ); // Iteration count from Matlab implementation
    EXPECT_TRUE( this->d_converged );

    // Compare against reference solution from Matlab
    std::vector<double> ref = {
        90.6666666666667,
       149.3333333333333,
       180.0000000000000,
       186.6666666666666,
       173.3333333333333,
       144.0000000000000,
       102.6666666666666,
        53.3333333333333};

    linalg_traits::test_vector<MV>(this->d_x,ref);

    // Solve again, should return without iterating
    this->solve();

    EXPECT_EQ( 0, this->d_iters );
    EXPECT_TRUE( this->d_converged );

    // Make sure solution didn't change
    linalg_traits::test_vector<MV>(this->d_x,ref);
}

//---------------------------------------------------------------------------//
//                        end of tstRichardson.cc
//---------------------------------------------------------------------------//
