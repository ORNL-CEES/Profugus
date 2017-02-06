//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/Utility_Functions_Tester.hh
 * \author Tom Evans
 * \date   Mon Dec 05 14:42:30 2016
 * \brief  Utility_Functions_Tester class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_test_Utility_Functions_Tester_hh
#define CudaUtils_test_Utility_Functions_Tester_hh

#include <vector>

#include "gtest/Gtest_Functions.hh"

//---------------------------------------------------------------------------//
// THREAD ID TESTS
//---------------------------------------------------------------------------//

class ThreadID_Test : public ::testing::Test
{
  protected:

    void run( const unsigned int (&nb)[2], const unsigned int (&nt)[3]);

    void test_1D_1D();
    void test_1D_2D();
    void test_1D_3D();
    void test_2D_1D();
    void test_2D_2D();
    void test_2D_3D();

    void test()
    {
        EXPECT_VEC_EQ(rtx, vtx);
        EXPECT_VEC_EQ(rty, vty);
        EXPECT_VEC_EQ(rtz, vtz);
        EXPECT_VEC_EQ(rbx, vbx);
        EXPECT_VEC_EQ(rby, vby);
    }

  protected:

    std::vector<int> rtx, rty, rtz, rbx, rby;
    std::vector<int> vtx, vty, vtz, vbx, vby;
};

#endif // CudaUtils_test_Utility_Functions_Tester_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Utility_Functions_Tester.hh
//---------------------------------------------------------------------------//
