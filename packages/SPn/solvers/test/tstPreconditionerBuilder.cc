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
    typedef profugus::PreconditionerBuilder<OP> Builder;
    typedef Builder::RCP_ParameterList          RCP_ParameterList;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Parallelism
        node  = profugus::node();
        nodes = profugus::nodes();

        // Build an Epetra map
        int global_size = 20;
        d_A = linalg_traits::build_matrix<Matrix>("laplacian",global_size);

        int my_size = d_A->NumMyRows();

        // Build lhs and rhs vectors
        d_x = linalg_traits::build_vector<MV>(global_size);
        d_y = linalg_traits::build_vector<MV>(global_size);
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_A->GRID(my_row);
            d_x->ReplaceGlobalValue(global_row,0,
                    static_cast<double>(20-global_row));
        }
    }

    void build_preconditioner( RCP_ParameterList db )
    {
        d_P = Builder::build_preconditioner(d_A,db);
    }

  protected:
    int node;
    int nodes;

    Teuchos::RCP<Epetra_CrsMatrix>       d_A;
    Teuchos::RCP<Epetra_Operator>        d_P;
    Teuchos::RCP<Epetra_MultiVector>     d_x;
    Teuchos::RCP<Epetra_MultiVector>     d_y;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST_F(PreconditionerBuilderTest, basic)
{
    RCP_ParameterList db = rcp(new Teuchos::ParameterList("test_db"));
    int ierr;
    std::vector<double> y_norm(1);

    // Default preconditioner should be valid
    std::cout << "Building default preconditioner" << std::endl;
    build_preconditioner(db);
    EXPECT_FALSE( d_P == Teuchos::null );
    d_y->PutScalar(0.0);
    ierr = d_P->Apply(*d_x,*d_y);
    EXPECT_FALSE(ierr);
    d_y->Norm2(&y_norm[0]);
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
    d_y->PutScalar(0.0);
    ierr = d_P->Apply(*d_x,*d_y);
    EXPECT_FALSE(ierr);
    d_y->Norm2(&y_norm[0]);
    EXPECT_TRUE(y_norm[0] > 1.0);

#ifdef USE_ML
    // Change preconditioner to "ML", should be valid
    std::cout << "Building ML prec" << std::endl;
    db->set("Preconditioner", std::string("ML"));
    build_preconditioner(db);
    EXPECT_FALSE( d_P == Teuchos::null );
    d_y->PutScalar(0.0);
    ierr = d_P->Apply(*d_x,*d_y);
    EXPECT_FALSE(ierr);
    d_y->Norm2(&y_norm[0]);
    EXPECT_TRUE(y_norm[0] > 1.0);
#endif
}

//---------------------------------------------------------------------------//
//                        end of tstPreconditionerBuilder.cc
//---------------------------------------------------------------------------//
