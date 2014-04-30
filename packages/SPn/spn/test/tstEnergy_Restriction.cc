//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstEnergy_Restriction.cc
 * \author Steven Hamilton
 * \date   Mon Apr 01 12:49:01 2013
 * \brief  Energy Grid Transfer test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <SPn/config.h>

#include <vector>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "../Energy_Restriction.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST(RestrictTest, Even)
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
    profugus::Energy_Restriction restrict0( vec0, vec1, steer );

    double tol=1.e-12;

    // Test restriction
    Epetra_Vector &fine_vec   = *(vec0(0));
    Epetra_Vector &coarse_vec = *(vec1(0));
    for( int i=0; i<fine_vec.MyLength(); ++i )
    {
        fine_vec[i] = static_cast<double>(2*i);
    }
    int error = restrict0.Apply(fine_vec,coarse_vec);
    for( int i=0; i<coarse_vec.MyLength(); ++i )
    {
        EXPECT_SOFTEQ(static_cast<double>(4*i+1),coarse_vec[i],tol);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Restriction.cc
//---------------------------------------------------------------------------//
