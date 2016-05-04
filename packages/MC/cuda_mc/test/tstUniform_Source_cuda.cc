//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstUniform_Source_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 15:22:15 2016
 * \brief  Test for Uniform_Source
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(UniformSource, coincident)
{
    Uniform_Source_Tester::test_source();
}

TEST(UniformSource, host_api)
{
    Uniform_Source_Tester::test_host_api();
}

//---------------------------------------------------------------------------//
//                 end of tstUniform_Source.cc
//---------------------------------------------------------------------------//
