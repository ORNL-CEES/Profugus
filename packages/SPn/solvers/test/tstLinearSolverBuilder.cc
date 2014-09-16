//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstLinearSolverBuilder.cc
 * \author Steven Hamilton
 * \brief  LinearSolverBuilder unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <SPn/config.h>

#include "../LinearSolverBuilder.hh"
#include "../Richardson.hh"
#include "../StratimikosSolver.hh"
#include "LinAlgTraits.hh"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

template <class T>
class SolverBuilderTest : public testing::Test
{
  protected:

    typedef typename T::MATRIX MATRIX;
    typedef profugus::LinearSolverBuilder<T> Builder;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
    }

    void build_solver( Teuchos::RCP<Teuchos::ParameterList> db )
    {
        d_solver = Builder::build_solver(db);
    }

  protected:

    Teuchos::RCP<profugus::LinearSolver<T> > d_solver;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(SolverBuilderTest, MyTypes);

TYPED_TEST(SolverBuilderTest, basic)
{
    Teuchos::RCP<Teuchos::ParameterList> db =
        Teuchos::rcp(new Teuchos::ParameterList("test_db"));

    // Default solver is Richardson
    this->build_solver(db);
    EXPECT_EQ("Profugus Richardson", this->d_solver->solver_label());
    Teuchos::RCP<profugus::Richardson<TypeParam> > rich =
        Teuchos::rcp_dynamic_cast<profugus::Richardson<TypeParam> >(this->d_solver);
    EXPECT_TRUE( rich != Teuchos::null );

    //
    // Profugus solver by specifying profugus_solver
    //

    db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));
    db->set("profugus_solver", std::string("Richardson"));
    this->build_solver(db);
    EXPECT_EQ("Profugus Richardson", this->d_solver->solver_label());
    rich = Teuchos::rcp_dynamic_cast<profugus::Richardson<TypeParam> >(this->d_solver);
    EXPECT_TRUE( rich != Teuchos::null );

    //
    // Profugus solver by specifying solver_type
    //

    db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));
    db->set("solver_type", std::string("Profugus"));
    this->build_solver(db);
    EXPECT_EQ("Profugus Richardson", this->d_solver->solver_label());
    rich = Teuchos::rcp_dynamic_cast<profugus::Richardson<TypeParam> >(this->d_solver);
    EXPECT_TRUE( rich != Teuchos::null );

    //
    // Stratimikos solver (default is AztecOO)
    //

    db = Teuchos::rcp(new Teuchos::ParameterList("test_db"));
    db->set("solver_type", std::string("Stratimikos"));
    this->build_solver(db);
    EXPECT_EQ("Stratimikos AztecOO", this->d_solver->solver_label());
    Teuchos::RCP<profugus::StratimikosSolver<TypeParam> > strat =
        Teuchos::rcp_dynamic_cast<profugus::StratimikosSolver<TypeParam> >(
            this->d_solver);
    EXPECT_TRUE( strat != Teuchos::null );
}

//---------------------------------------------------------------------------//
//                        end of tstLinearSolverBuilder.cc
//---------------------------------------------------------------------------//
