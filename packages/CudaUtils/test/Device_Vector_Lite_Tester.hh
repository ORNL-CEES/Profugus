//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/Device_Vector_Lite_Tester.hh
 * \author Steven Hamilton
 * \date   Tue May 03 12:09:58 2016
 * \brief  Device_Vector_Lite_Tester class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_test_Device_Vector_Lite_Tester_hh
#define CudaUtils_test_Device_Vector_Lite_Tester_hh

//===========================================================================//
/*!
 * \class Device_Vector_Lite_Tester
 * \brief Test Device_Vector_Lite on-device functions.
 */
//===========================================================================//

class Device_Vector_Lite_Tester
{
  public:

    static void test_host();
    static void test_device();
    static void test_host_copy();
};

//---------------------------------------------------------------------------//
#endif // CudaUtils_test_Device_Vector_Lite_Tester_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Device_Vector_Lite_Tester.hh
//---------------------------------------------------------------------------//
