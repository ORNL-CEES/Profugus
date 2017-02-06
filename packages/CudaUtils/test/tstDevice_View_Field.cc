//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/tstDevice_View_Field.cc
 * \author Steven Hamilton
 * \date   Fri Nov 18 11:55:39 2016
 * \brief  Tests for class Device_View_Field.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Device_View_Field_Tester.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(DeviceViewFieldTest, views)
{
    Device_View_Field_Tester::test_views();
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/tstDevice_View_Field.cc
//---------------------------------------------------------------------------//
