//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstPartitioner.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 12 09:55:23 2014
 * \brief  Partitioner unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Partitioner.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class Partitioner_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::Partitioner Partitioner;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> Data that get re-initialized between tests
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Partitioner_Test, stub)
{
}

//---------------------------------------------------------------------------//
//                 end of tstPartitioner.cc
//---------------------------------------------------------------------------//
