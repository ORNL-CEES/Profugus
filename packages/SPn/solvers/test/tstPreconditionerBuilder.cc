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

using Teuchos::RCP;
using Teuchos::rcp;

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class PreconditionerBuilderTest : public ::testing::Test
{
  protected:

    typedef profugus::PreconditionerBuilder Builder;
    typedef Builder::RCP_ParameterList      RCP_ParameterList;
    typedef Epetra_MultiVector              MV;
    typedef Epetra_Operator                 OP;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // Parallelism
        node  = profugus::node();
        nodes = profugus::nodes();

        // Build Epetra communicator
#ifdef COMM_MPI
        Epetra_MpiComm comm(profugus::communicator);
#else
        Epetra_SerialComm comm;
#endif

        // Build an Epetra map
        int global_size = 20;
        d_map = Teuchos::rcp( new Epetra_Map( global_size, 0, comm ) );

        int my_size = global_size / nodes;

        // Build CrsMatrix
        d_A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,*d_map,3) );
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_map->GID(my_row);
            if( global_row == 0 )
            {
                std::vector<int> ids(2);
                ids[0] = 0;
                ids[1] = 1;
                std::vector<double> vals(2);
                vals[0] = 1.0;
                vals[1] = 1.0;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else if( global_row == global_size-1 )
            {
                std::vector<int> ids(2);
                ids[0] = 18;
                ids[1] = 19;
                std::vector<double> vals(2);
                vals[0] = -1.0;
                vals[1] = 20.0;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else
            {
                std::vector<int> ids(3);
                ids[0] = global_row-1;
                ids[1] = global_row;
                ids[2] = global_row+1;
                std::vector<double> vals(3);
                vals[0] = -1.0;
                vals[1] = static_cast<double>(global_row+1);
                vals[2] =  1.0;
                d_A->InsertGlobalValues(global_row,3,&vals[0],&ids[0]);
            }
        }
        d_A->FillComplete();

        // Build lhs and rhs vectors
        d_x = Teuchos::rcp( new Epetra_MultiVector(*d_map,1) );
        d_y = Teuchos::rcp( new Epetra_MultiVector(*d_map,1) );
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_map->GID(my_row);
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

    Teuchos::RCP<Epetra_Map>             d_map;
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
