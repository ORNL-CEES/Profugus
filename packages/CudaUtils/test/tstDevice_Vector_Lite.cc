//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/tstDevice_Vector_Lite.cc
 * \author Steven Hamilton
 * \date   Tue May 03 12:04:39 2016
 * \brief  Device_Vector_Lite class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <type_traits>


#include "Device_Vector_Lite_Tester.hh"
#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(DeviceVectorLite, host)
{
    Device_Vector_Lite_Tester::test_host();
}

TEST(DeviceVectorLite, device)
{
    Device_Vector_Lite_Tester::test_device();
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/tstDevice_Vector_Lite.cc
//---------------------------------------------------------------------------//
