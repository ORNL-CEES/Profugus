//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstDomain_Transporter_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Domain_Transporter
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Domain_Transporter_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(DomainTransporter, three_group)
{
    Domain_Transporter_Tester::test_transport(3);
}

TEST(DomainTransporter, five_group)
{
    Domain_Transporter_Tester::test_transport(5);
}

//---------------------------------------------------------------------------//
//                 end of tstDomain_Transporter_cuda.cc
//---------------------------------------------------------------------------//
