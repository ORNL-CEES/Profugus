//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFission_Matrix_Acceleration.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 12 14:59:19 2014
 * \brief  Test for Fission_Matrix_Acceleration
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Fission_Matrix_Acceleration.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class FM_AccelerationTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Fission_Matrix_Acceleration Acceleration;

  protected:
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(FM_AccelerationTest, simple)
{
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Matrix_Acceleration.cc
//---------------------------------------------------------------------------//
