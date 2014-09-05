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

template <class T>
class PreconditionerBuilderTest : public ::testing::Test
{
  protected:

    typedef typename linalg_traits::traits_types<T>::MV     MV;
    typedef typename linalg_traits::traits_types<T>::OP     OP;
    typedef typename linalg_traits::traits_types<T>::Matrix Matrix;

    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef profugus::PreconditionerBuilder<OP>   Builder;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Build an Epetra map
        d_N = 20;
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",d_N);

        // Build lhs and rhs vectors
        d_x = linalg_traits::build_vector<MV>(d_N);
        d_y = linalg_traits::build_vector<MV>(d_N);
        std::vector<double> vals(d_N);
        for( int i=0; i<d_N; ++i )
            vals[i] = static_cast<double>(d_N-i);
        linalg_traits::fill_vector<MV>(d_x,vals);
    }

    void build_preconditioner( RCP<Teuchos::ParameterList> db )
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
typedef ::testing::Types<Epetra_MultiVector,Tpetra_MultiVector> MyTypes;
TYPED_TEST_CASE(PreconditionerBuilderTest, MyTypes);

TYPED_TEST(PreconditionerBuilderTest, basic)
{
    typedef typename TestFixture::MV  MV;
    typedef typename TestFixture::OPT OPT;
    typedef typename TestFixture::MVT MVT;

    RCP<Teuchos::ParameterList> db =
        rcp(new Teuchos::ParameterList("test_db"));
    std::vector<double> y_norm(1);

    // Determine if this is epetra or tpetra
    Teuchos::RCP<Epetra_MultiVector> x_ep =
        Teuchos::rcp_dynamic_cast<Epetra_MultiVector>(this->d_x);
    bool epetra = (x_ep != Teuchos::null);

    // Default preconditioner should be valid
    std::cout << "Building default preconditioner" << std::endl;
    this->build_preconditioner(db);
    EXPECT_FALSE( this->d_P == Teuchos::null );
    std::vector<double> zero(this->d_N,0.0);
    linalg_traits::fill_vector<MV>(this->d_y,zero);
    OPT::Apply(*this->d_P,*this->d_x,*this->d_y);
    MVT::MvNorm(*this->d_y,y_norm);
    EXPECT_TRUE(y_norm[0] > 1.0);

    // Set preconditioner to "None", should return null RCP
    std::cout << "Building null prec" << std::endl;
    db->set("Preconditioner", std::string("None"));
    this->build_preconditioner(db);
    EXPECT_TRUE( this->d_P == Teuchos::null );

    if( epetra )
    {
        // Change preconditioner to "Ifpack", should be valid for epetra
        std::cout << "Building Ifpack prec" << std::endl;
        db->set("Preconditioner", std::string("Ifpack"));
        this->build_preconditioner(db);
        EXPECT_FALSE( this->d_P == Teuchos::null );
        linalg_traits::fill_vector<MV>(this->d_y,zero);
        OPT::Apply(*this->d_P,*this->d_x,*this->d_y);
        MVT::MvNorm(*this->d_y,y_norm);
        EXPECT_TRUE(y_norm[0] > 1.0);
    }

#ifdef USE_ML
    if( epetra )
    {
        // Change preconditioner to "ML", should be valid
        std::cout << "Building ML prec" << std::endl;
        db->set("Preconditioner", std::string("ML"));
        this->build_preconditioner(db);
        EXPECT_FALSE( this->d_P == Teuchos::null );
        linalg_traits::fill_vector<MV>(this->d_y,zero);
        OPT::Apply(*this->d_P,*this->d_x,*this->d_y);
        MVT::MvNorm(*this->d_y,y_norm);
        EXPECT_TRUE(y_norm[0] > 1.0);
    }
#endif

    if( !epetra )
    {
        // Change preconditioner to "Ifpack2", should be valid for tpetra
        std::cout << "Building Ifpack2 prec" << std::endl;
        db->set("Preconditioner", std::string("Ifpack2"));
        this->build_preconditioner(db);
        EXPECT_FALSE( this->d_P == Teuchos::null );
        linalg_traits::fill_vector<MV>(this->d_y,zero);
        OPT::Apply(*this->d_P,*this->d_x,*this->d_y);
        MVT::MvNorm(*this->d_y,y_norm);
        EXPECT_TRUE(y_norm[0] > 1.0);
    }

    if( !epetra )
    {
        // Change preconditioner to "MueLu", should be valid for tpetra
        std::cout << "Building MueLu prec" << std::endl;
        db->set("Preconditioner", std::string("MueLu"));
        this->build_preconditioner(db);
        EXPECT_FALSE( this->d_P == Teuchos::null );
        linalg_traits::fill_vector<MV>(this->d_y,zero);
        OPT::Apply(*this->d_P,*this->d_x,*this->d_y);
        MVT::MvNorm(*this->d_y,y_norm);
        EXPECT_TRUE(y_norm[0] > 1.0);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstPreconditionerBuilder.cc
//---------------------------------------------------------------------------//
