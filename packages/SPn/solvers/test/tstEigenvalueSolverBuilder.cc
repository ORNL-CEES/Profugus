//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/test/tstEigenvalueSolverBuilder.cc
 * \author Steven Hamilton
 * \brief  EigenvalueSolverBuilder unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <SPn/config.h>

#include "../EigenvalueSolverBuilder.hh"
#include "../Arnoldi.hh"
#include "../Davidson_Eigensolver.hh"
#include "../PowerIteration.hh"
#include "../RayleighQuotient.hh"
#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

template <class T>
class EigenvalueSolverBuilderTest : public ::testing::Test
{
  protected:

    typedef typename T::MV     MV;
    typedef typename T::OP     OP;
    typedef typename T::MATRIX MATRIX;

    typedef profugus::EigenvalueSolverBuilder<T> Builder;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        int num_global = 4;
        d_A = linalg_traits::build_matrix<T>("laplacian",num_global);
        d_B = linalg_traits::build_matrix<T>("diagonal",num_global);
    }

  protected:

    Teuchos::RCP<profugus::EigenvalueSolver<T> > d_solver;
    Teuchos::RCP<MATRIX> d_A;
    Teuchos::RCP<MATRIX> d_B;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(EigenvalueSolverBuilderTest, MyTypes);

TYPED_TEST(EigenvalueSolverBuilderTest, basic)
{
    typedef profugus::EigenvalueSolverBuilder<TypeParam> Builder;

    Teuchos::RCP<Teuchos::ParameterList> db =
        Teuchos::rcp(new Teuchos::ParameterList("test_db"));

    // Default standard eigenvalue solver is Arnoldi
    this->d_solver = Builder::build_solver(db,this->d_A);
    EXPECT_EQ("Arnoldi",this->d_solver->solver_label());
    Teuchos::RCP<profugus::Arnoldi<TypeParam> > arnoldi =
        Teuchos::rcp_dynamic_cast<profugus::Arnoldi<TypeParam> >(this->d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    // Default generalized eigenvalue solver is Arnoldi (for now)
    this->d_solver = Builder::build_solver(db,this->d_A,this->d_B);
    EXPECT_EQ("Arnoldi",this->d_solver->solver_label());
    arnoldi = Teuchos::rcp_dynamic_cast<profugus::Arnoldi<TypeParam> >(this->d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    // Make sure "Arnoldi" keyword is recognized by both functions
    db->set("eigensolver",std::string("Arnoldi"));

    this->d_solver = Builder::build_solver(db,this->d_A);
    EXPECT_EQ("Arnoldi",this->d_solver->solver_label());
    arnoldi = Teuchos::rcp_dynamic_cast<profugus::Arnoldi<TypeParam> >(this->d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    this->d_solver = Builder::build_solver(db,this->d_A,this->d_B);
    EXPECT_EQ("Arnoldi",this->d_solver->solver_label());
    arnoldi = Teuchos::rcp_dynamic_cast<profugus::Arnoldi<TypeParam> >(this->d_solver);
    EXPECT_TRUE( arnoldi != Teuchos::null );

    // Power iteration for both standard and generalized problems
    db->set("eigensolver",std::string("Power"));

    this->d_solver = Builder::build_solver(db,this->d_A);
    EXPECT_EQ("Power Iteration",this->d_solver->solver_label());
    Teuchos::RCP<profugus::PowerIteration<TypeParam> > power =
        Teuchos::rcp_dynamic_cast<profugus::PowerIteration<TypeParam> >(this->d_solver);
    EXPECT_TRUE( power != Teuchos::null );

    this->d_solver = Builder::build_solver(db,this->d_A,this->d_B);
    EXPECT_EQ("Power Iteration",this->d_solver->solver_label());
    power = Teuchos::rcp_dynamic_cast<profugus::PowerIteration<TypeParam> >(this->d_solver);
    EXPECT_TRUE( power != Teuchos::null );

    // Rayleigh quotient iteration for generalized problem
    db->set("eigensolver",std::string("RQI"));

    this->d_solver = Builder::build_solver(db,this->d_A,this->d_B);
    EXPECT_EQ("Rayleigh Quotient",this->d_solver->solver_label());
    Teuchos::RCP<profugus::RayleighQuotient<TypeParam> > rqi =
        Teuchos::rcp_dynamic_cast<profugus::RayleighQuotient<TypeParam> >(this->d_solver);
    EXPECT_TRUE( rqi != Teuchos::null );

    // Davidson for generalized problem
    db->set("eigensolver",std::string("Davidson"));

    this->d_solver = Builder::build_solver(db,this->d_A,this->d_B);
    EXPECT_EQ("Davidson",this->d_solver->solver_label());
    Teuchos::RCP<profugus::Davidson_Eigensolver<TypeParam> > davidson =
        Teuchos::rcp_dynamic_cast<profugus::Davidson_Eigensolver<TypeParam> >(this->d_solver);
    EXPECT_TRUE( davidson!= Teuchos::null );
}

//---------------------------------------------------------------------------//
//                        end of tstEigenvalueSolverBuilder.cc
//---------------------------------------------------------------------------//
