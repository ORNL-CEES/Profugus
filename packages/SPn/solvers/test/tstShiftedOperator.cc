//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstShiftedOperator.cc
 * \author Steven Hamilton
 * \brief  ShiftedOperator unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../ShiftedOperator.hh"

#include <SPn/config.h>

#include "comm/global.hh"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

//---------------------------------------------------------------------------//
// Test fixture base class
//---------------------------------------------------------------------------//

class ShiftedOperatorTest : public testing::Test
{
  protected:
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
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
        d_B = Teuchos::rcp( new Epetra_CrsMatrix(Copy,*d_map,1) );
        for( int my_row=0; my_row<my_size; ++my_row )
        {
            int global_row = d_map->GID(my_row);
            if( global_row == 0 )
            {
                std::vector<int> ids(2);
                ids[0] = 0;
                ids[1] = 1;
                std::vector<double> vals(2);
                vals[0] =  2.0;
                vals[1] = -1.0;
                d_A->InsertGlobalValues(global_row,2,&vals[0],&ids[0]);
            }
            else if( global_row == global_size-1 )
            {
                std::vector<int> ids(2);
                ids[0] = 18;
                ids[1] = 19;
                std::vector<double> vals(2);
                vals[0] = -1.0;
                vals[1] =  2.0;
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
                vals[1] =  2.0;
                vals[2] = -1.0;
                d_A->InsertGlobalValues(global_row,3,&vals[0],&ids[0]);
            }
            std::vector<int>    inds(1);
            std::vector<double> vals(1);
            inds[0] = global_row;
            vals[0] = static_cast<double>(global_row+1);
            d_B->InsertGlobalValues(global_row,1,&vals[0],&inds[0]);
        }
        d_A->FillComplete();
        d_B->FillComplete();

        // Build eigenvector
        d_x = Teuchos::rcp( new Epetra_Vector(*d_map) );

        // Build solver
        d_operator = Teuchos::rcp(new profugus::ShiftedOperator<MV,OP>());
        Check(!d_operator.is_null());
        d_operator->set_operator(d_A);
        d_operator->set_rhs_operator(d_B);
    }

  protected:
    int node;
    int nodes;

    Teuchos::RCP<Epetra_Map>                        d_map;
    Teuchos::RCP<Epetra_CrsMatrix>                  d_A;
    Teuchos::RCP<Epetra_CrsMatrix>                  d_B;
    Teuchos::RCP<Epetra_Vector>                     d_x;
    Teuchos::RCP<profugus::ShiftedOperator<MV,OP> > d_operator;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ShiftedOperatorTest, basic)
{
    //Clone vector for holding solutions
    Teuchos::RCP<Epetra_Vector> y(new Epetra_Vector(*d_map));

    // Unshifted operator, same as multiplying by A
    d_operator->set_shift(0.0);
    d_x->PutScalar(1.0);
    d_operator->Apply(*d_x,*y);

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_map->GID(i);
        if( global_row == 0 || global_row == 19 )
        {
            EXPECT_SOFTEQ( 1.0, (*y)[i], 1e-14 );
        }
        else
        {
            EXPECT_SOFTEQ( 0.0, (*y)[i], 1e-14 );
        }
    }

    // New vector
    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_map->GID(i);
        (*d_x)[i] = static_cast<double>(global_row+1);
    }

    d_operator->Apply(*d_x,*y);
    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_map->GID(i);
        if( global_row == 19 )
        {
            EXPECT_SOFTEQ( 21.0, (*y)[i], 1e-14 );
        }
        else
        {
            EXPECT_SOFTEQ( 0.0, (*y)[i], 1e-14 );
        }
    }

    // Now set a shift
    d_operator->set_shift(0.5);
    d_x->PutScalar(1.0);
    d_operator->Apply(*d_x,*y);

    // Matlab computed reference
    double ref[] = {
        0.5000,  -1.0000,  -1.5000,  -2.0000,  -2.5000,  -3.0000,  -3.5000,
       -4.0000,  -4.5000,  -5.0000,  -5.5000,  -6.0000,  -6.5000,  -7.0000,
       -7.5000,  -8.0000,  -8.5000,  -9.0000,  -9.5000,  -9.0000 };

    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_map->GID(i);
        EXPECT_SOFTEQ( ref[global_row], (*y)[i], 1e-14 );
    }


    // Different vector
    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_map->GID(i);
        (*d_x)[i] = static_cast<double>(global_row+1);
    }
    d_operator->Apply(*d_x,*y);


    // Matlab computed reference
    double ref2[] = {
        -0.5000,   -2.0000,   -4.5000,   -8.0000,  -12.5000, -18.0000,
       -24.5000,  -32.0000,  -40.5000,  -50.0000,  -60.5000, -72.0000,
       -84.5000,  -98.0000, -112.5000, -128.0000, -144.5000, -162.0000,
      -180.5000, -179.0000 };
    for( int i=0; i<y->MyLength(); ++i )
    {
        int global_row = d_map->GID(i);
        EXPECT_SOFTEQ( ref2[global_row], (*y)[i], 1e-14 );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstShiftedOperator.cc
//---------------------------------------------------------------------------//
