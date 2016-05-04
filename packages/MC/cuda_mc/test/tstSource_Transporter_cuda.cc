//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstSource_Transporter_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Source_Transporter
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source_Transporter_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Source_Transporter_cudaTest, onegroup)
{
    Source_Transporter_Tester::test_transport(1);
}

TEST(Source_Transporter_cudaTest, threegroup)
{
    Source_Transporter_Tester::test_transport(3);
}

TEST(Source_Transporter_cudaTest, fivegroup)
{
    Source_Transporter_Tester::test_transport(5);
}

//---------------------------------------------------------------------------//
//                 end of tstSource_Transporter_cuda.cc
//---------------------------------------------------------------------------//
