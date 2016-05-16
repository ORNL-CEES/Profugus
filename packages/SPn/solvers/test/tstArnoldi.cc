//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/test/tstArnoldi.cc
 * \author Steven Hamilton
 * \date   Mon Aug 22 12:22:41 2011
 * \brief  Test Arnoldi class.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <random>

#include <SPn/config.h>

#include "../Arnoldi.hh"
#include "../InverseOperator.hh"
#include "LinAlgTraits.hh"

#include "Teuchos_RCP.hpp"

//---------------------------------------------------------------------------//
// Test Fixture
//---------------------------------------------------------------------------//

template <class T>
class Arnoldi_Test : public testing::Test
{
  protected:
    void SetUp()
    {
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(Arnoldi_Test, MyTypes);

TYPED_TEST(Arnoldi_Test, Eigensolver)
{
    typedef typename TypeParam::MV     MV;
    typedef typename TypeParam::MATRIX MATRIX;

    typedef profugus::Arnoldi<TypeParam>         Arnoldi;
    typedef profugus::InverseOperator<TypeParam> InverseOperator;

    using Teuchos::RCP;
    RCP<Teuchos::ParameterList> db(new Teuchos::ParameterList("test"));
    RCP<Teuchos::ParameterList> db2(new Teuchos::ParameterList("test"));

    // Matrix size
    int N = 8;

    RCP<MATRIX> A = linalg_traits::build_matrix<TypeParam>("laplacian",N);
    RCP<MATRIX> B = linalg_traits::build_matrix<TypeParam>("diagonal",N);
    RCP<MV> evec = linalg_traits::build_vector<TypeParam>(N);
    std::random_device rd;
    std::vector<double> init(N);
    for( int i=0; i<N; ++i )
        init[i] = rd();
    linalg_traits::fill_vector<TypeParam>(evec,init);

    // Test with matrix A
    double e_val;

    // database settings
    {
        // make anasazi settings
        Teuchos::ParameterList adb("anasazi");
        adb.set("Convergence Tolerance", 1.0e-10);
        adb.set("Maximum Restarts", 50);
        db->set("Anasazi", adb);
    }

    // make solver
    Arnoldi solver( db );
    solver.set_operator(A);
    // This is actually not a great test problem because the dominant
    // eigenvector is exactly orthogonal to a constant vector.  We have
    // to initialize with a random vector or else we'll converge to the
    // 2nd eigenmode.
    solver.solve( e_val, evec );

    std::cout.precision(15);
    std::cout << "Eig(A) = " << e_val << endl;

    std::vector<double> eref =
      { 0.161229841765317, -0.303012985114695, 0.408248290463863,
       -0.464242826880012, 0.464242826880013, -0.408248290463863,
        0.303012985114696, -0.161229841765317};

    // Test against Matlab computed values
    double tol = 1.0e-6;
    EXPECT_SOFTEQ(e_val, 3.879385241571816, tol);

    linalg_traits::set_sign<TypeParam>(evec);
    linalg_traits::test_vector<TypeParam>(evec,eref);

    // Create operator for A^{-1}B
    {
        db2->set("max_itr", 50);
        db2->set("tolerance", 1e-12);
        db2->set("solver_type", std::string("stratimikos"));
    }

    RCP<InverseOperator> AinvB(Teuchos::rcp( new InverseOperator(db2)));
    AinvB->set_operator(A);
    AinvB->set_rhs_operator(B);

    // Now an eigensolver for A^{-1}B
    Arnoldi gensolver( db );
    gensolver.set_operator( AinvB );

    // Solve
    evec = linalg_traits::build_vector<TypeParam>(N);
    gensolver.solve( e_val, evec );

    std::cout << "Eig(AinvB) = " << e_val << endl;

    std::vector<double> eref2 =
        {0.125730096111867, 0.248228875407670, 0.357968479865590,
         0.440108258878819, 0.477004161945638, 0.452604018262593,
         0.358411182752022, 0.199739110973060};

    // Test against Matlab computed values
    EXPECT_SOFTEQ(e_val, 38.909863460868493, tol);

    linalg_traits::set_sign<TypeParam>(evec);
    linalg_traits::test_vector<TypeParam>(evec,eref2);
}

//---------------------------------------------------------------------------//
//                        end of tstArnoldi.cc
//---------------------------------------------------------------------------//
