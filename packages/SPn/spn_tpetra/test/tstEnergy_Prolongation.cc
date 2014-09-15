//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/test/tstEnergy_Prolongation.cc
 * \author Steven Hamilton
 * \date   Mon Apr 01 12:49:01 2013
 * \brief  Energy Grid Transfer test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <SPn/config.h>

#include <vector>

#include "Teuchos_DefaultComm.hpp"

#include "solvers/LinAlgTypedefs.hh"
#include "../Energy_Prolongation.hh"

typedef typename LinAlgTypedefs<TPETRA>::MV     MV;
typedef typename LinAlgTypedefs<TPETRA>::OP     OP;
typedef typename LinAlgTypedefs<TPETRA>::VECTOR VECTOR;
typedef typename LinAlgTypedefs<TPETRA>::MAP    MAP;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST(ProlongTest, Even)
{
    int Nv = 50;
    int Ng = 8;

    // Get Teuchos communicator
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    // Create Epetra maps
    Teuchos::RCP<MAP> map0( new MAP(
                Teuchos::OrdinalTraits<int>::invalid(),Nv*Ng  ,0,comm) );
    Teuchos::RCP<MAP> map1( new MAP(
                Teuchos::OrdinalTraits<int>::invalid(),Nv*Ng/2,0,comm) );

    // Create Epetra vectors
    Teuchos::RCP<MV> vec0( new MV(map0,1) );
    Teuchos::RCP<MV> vec1( new MV(map1,1) );

    std::vector<int> steer(4,2);
    profugus::tpetra::Energy_Prolongation prolong0( vec1, vec0, steer );

    double tol=1.e-12;

    // Test prolongation
    Teuchos::ArrayRCP<double> fine_data   = vec0->getDataNonConst(0);
    Teuchos::ArrayRCP<double> coarse_data = vec1->getDataNonConst(0);
    for( int i=0; i<vec1->getLocalLength(); ++i )
    {
        coarse_data[i] = static_cast<double>(i);
    }
    prolong0.apply(*vec1,*vec0);
    for( int i=0; i<vec0->getLocalLength(); ++i )
    {
        EXPECT_SOFTEQ(static_cast<double>(i/2),fine_data[i],tol);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Prolongation.cc
//---------------------------------------------------------------------------//
