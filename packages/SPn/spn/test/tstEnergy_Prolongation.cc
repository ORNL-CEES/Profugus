//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstEnergy_Prolongation.cc
 * \author Steven Hamilton
 * \date   Mon Apr 01 12:49:01 2013
 * \brief  Energy Grid Transfer test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>

#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "solvers/LinAlgTypedefs.hh"
#include "../Energy_Prolongation.hh"
#include "../MatrixTraits.hh"
#include "../VectorTraits.hh"


//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class T>
class ProlongTest : public testing::Test
{
  protected:
      void SetUp(){};
};

using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(ProlongTest, MyTypes);

TYPED_TEST(ProlongTest, Even)
{
    typedef typename TypeParam::MAP    Map_t;
    typedef typename TypeParam::OP     OP;
    typedef typename TypeParam::MV     MV;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

    int Nv = 50;
    int Ng = 8;

    int nodes = profugus::nodes();

    // Create maps
    Teuchos::RCP<Map_t> map0 =
        profugus::MatrixTraits<TypeParam>::build_map(Nv*Ng,Nv*Ng*nodes);
    Teuchos::RCP<Map_t> map1 =
        profugus::MatrixTraits<TypeParam>::build_map(Nv*Ng/2,Nv*Ng*nodes/2);

    // Create Epetra vectors
    Teuchos::RCP<MV> vec0 =
        profugus::VectorTraits<TypeParam>::build_vector(map0);
    Teuchos::RCP<MV> vec1 =
        profugus::VectorTraits<TypeParam>::build_vector(map1);

    std::vector<int> steer(4,2);
    profugus::Energy_Prolongation<TypeParam> prolong0( map1, map0, steer );

    double tol=1.e-12;

    // Test prolongation
    Teuchos::ArrayView<double> fine_data =
        profugus::VectorTraits<TypeParam>::get_data_nonconst(vec0,0);
    Teuchos::ArrayView<double> coarse_data =
        profugus::VectorTraits<TypeParam>::get_data_nonconst(vec1,0);

    for( int i=0; i<coarse_data.size(); ++i )
    {
        coarse_data[i] = static_cast<double>(i);
    }

    OPT::Apply(prolong0,*vec1,*vec0);

    for( int i=0; i<fine_data.size(); ++i )
    {
        EXPECT_SOFTEQ(static_cast<double>(i/2),fine_data[i],tol);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Prolongation.cc
//---------------------------------------------------------------------------//
