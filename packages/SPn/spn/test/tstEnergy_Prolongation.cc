//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstEnergy_Prolongation.cc
 * \author Steven Hamilton
 * \date   Mon Apr 01 12:49:01 2013
 * \brief  Energy Grid Transfer test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <SPn/config.h>

#include <vector>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "../Energy_Prolongation.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST(ProlongTest, Even)
{
    int Nv = 50;
    int Ng = 8;

    // Define Epetra communicator
#ifdef COMM_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    // Create Epetra maps
    Epetra_Map map0(-1,Nv*Ng  ,0,comm);
    Epetra_Map map1(-1,Nv*Ng/2,0,comm);

    // Create Epetra vectors
    Epetra_MultiVector vec0(map0,1);
    Epetra_MultiVector vec1(map1,1);

    std::vector<int> steer(4,2);
    profugus::Energy_Prolongation prolong0( vec1, vec0, steer );

    double tol=1.e-12;

    // Test prolongation
    Epetra_Vector &fine_vec   = *(vec0(0));
    Epetra_Vector &coarse_vec = *(vec1(0));
    for( int i=0; i<coarse_vec.MyLength(); ++i )
    {
        coarse_vec[i] = static_cast<double>(i);
    }
    int error = prolong0.Apply(coarse_vec,fine_vec);
    for( int i=0; i<fine_vec.MyLength(); ++i )
    {
        EXPECT_SOFTEQ(static_cast<double>(i/2),fine_vec[i],tol);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Prolongation.cc
//---------------------------------------------------------------------------//
