//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstPreconditionerBuilder.cc
 * \author Steven Hamilton
 * \brief  PreconditionerBuilder unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <string>
#include <vector>

#include <SPn/config.h>

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "../PreconditionerBuilder.hh"

#include "LinAlgTraits.hh"

using Teuchos::RCP;
using Teuchos::rcp;

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class PreconditionerBuilderTest : public ::testing::Test
{
  protected:

    typedef Epetra_MultiVector                  MV;
    typedef Epetra_Operator                     OP;
    typedef Epetra_CrsMatrix                    Matrix;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;
    typedef Anasazi::MultiVecTraits<double,MV>  MVT;
    typedef profugus::PreconditionerBuilder<OP> Builder;
    typedef Builder::RCP_ParameterList          RCP_ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Build an Epetra map
        d_N = 20;
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",d_N);

        int my_size = d_A->NumMyRows();

        // Build lhs and rhs vectors
        d_x = linalg_traits::build_vector<MV>(d_N);
        d_y = linalg_traits::build_vector<MV>(d_N);
        std::vector<double> vals(d_N);
        for( int i=0; i<d_N; ++i )
            vals[i] = static_cast<double>(d_N-i);
        linalg_traits::fill_vector<MV>(d_x,vals);
    }

    void build_preconditioner( RCP_ParameterList db )
    {
        d_P = Builder::build_preconditioner(d_A,db);
    }

  protected:

    int d_N;

    Teuchos::RCP<Matrix> d_A;
    Teuchos::RCP<OP>     d_P;
    Teuchos::RCP<MV>     d_x;
    Teuchos::RCP<MV>     d_y;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST_F(PreconditionerBuilderTest, basic)
{
    RCP_ParameterList db = rcp(new Teuchos::ParameterList("test_db"));
    std::vector<double> y_norm(1);

    // Default preconditioner should be valid
    std::cout << "Building default preconditioner" << std::endl;
    build_preconditioner(db);
    EXPECT_FALSE( d_P == Teuchos::null );
    std::vector<double> zero(d_N,0.0);
    linalg_traits::fill_vector<MV>(d_y,zero);
    OPT::Apply(*d_P,*d_x,*d_y);
    MVT::MvNorm(*d_y,y_norm);
    EXPECT_TRUE(y_norm[0] > 1.0);

    // Set preconditioner to "None", should return null RCP
    std::cout << "Building null prec" << std::endl;
    db->set("Preconditioner", std::string("None"));
    build_preconditioner(db);
    EXPECT_TRUE( d_P == Teuchos::null );

    // Change preconditioner to "Ifpack", should be valid
    std::cout << "Building Ifpack prec" << std::endl;
    db->set("Preconditioner", std::string("Ifpack"));
    build_preconditioner(db);
    EXPECT_FALSE( d_P == Teuchos::null );
    linalg_traits::fill_vector<MV>(d_y,zero);
    OPT::Apply(*d_P,*d_x,*d_y);
    MVT::MvNorm(*d_y,y_norm);
    EXPECT_TRUE(y_norm[0] > 1.0);

#ifdef USE_ML
    // Change preconditioner to "ML", should be valid
    std::cout << "Building ML prec" << std::endl;
    db->set("Preconditioner", std::string("ML"));
    build_preconditioner(db);
    EXPECT_FALSE( d_P == Teuchos::null );
    linalg_traits::fill_vector<MV>(d_y,zero);
    OPT::Apply(*d_P,*d_x,*d_y);
    MVT::MvNorm(*d_y,y_norm);
    EXPECT_TRUE(y_norm[0] > 1.0);
#endif
}

//---------------------------------------------------------------------------//
//                        end of tstPreconditionerBuilder.cc
//---------------------------------------------------------------------------//
