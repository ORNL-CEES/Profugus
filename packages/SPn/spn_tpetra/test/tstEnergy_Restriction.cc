//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/test/tstEnergy_Restriction.cc
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
#include "../Energy_Restriction.hh"

using profugus::TpetraTypes;
typedef typename TpetraTypes::MV     MV;
typedef typename TpetraTypes::OP     OP;
typedef typename TpetraTypes::VECTOR VECTOR;
typedef typename TpetraTypes::MAP    MAP;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST(RestrictTest, Even)
{
    int Nv = 50;
    int Ng = 8;

    // Get Teuchos communicator
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    // Create Epetra maps
    Teuchos::RCP<MAP> map0(new MAP(
                Teuchos::OrdinalTraits<int>::invalid(), Nv*Ng  ,0,comm) );
    Teuchos::RCP<MAP> map1(new MAP(
                Teuchos::OrdinalTraits<int>::invalid(), Nv*Ng/2,0,comm) );

    // Create Epetra vectors
    Teuchos::RCP<MV> vec0( new MV(map0,1) );
    Teuchos::RCP<MV> vec1( new MV(map1,1) );

    std::vector<int> steer(4,2);
    profugus::tpetra::Energy_Restriction restrict0( vec0, vec1, steer );

    double tol=1.e-12;

    // Test restriction
    Teuchos::ArrayRCP<double> fine_data   = vec0->getDataNonConst(0);
    Teuchos::ArrayRCP<double> coarse_data = vec1->getDataNonConst(0);
    for( int i=0; i<vec0->getLocalLength(); ++i )
    {
        fine_data[i] = static_cast<double>(2*i);
    }
    restrict0.apply(*vec0,*vec1);
    for( int i=0; i<vec1->getLocalLength(); ++i )
    {
        EXPECT_SOFTEQ(static_cast<double>(4*i+1),coarse_data[i],tol);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Restriction.cc
//---------------------------------------------------------------------------//
