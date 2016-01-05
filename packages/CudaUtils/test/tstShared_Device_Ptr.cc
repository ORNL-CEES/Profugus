//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstShared_Device_Ptr.hh
 * \author Stuart Slattery
 * \date   Thu Dec 17 11:43:04 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "Shared_Device_Ptr_Tester.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class SharedDevicePtrTest : public ::testing::Test
{

  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// Test pointer constructor.
TEST_F(SharedDevicePtrTest, pointer_test)
{
    // Test data.
    double v1 = 3.432l;
    int v2 = 4;

    // Make a foo.
    std::shared_ptr<Foo> foo = std::make_shared<Foo>( v1, v2 );

    // Make some bars.
    Bar b1( foo );
    Bar b2 = b1;
    Bar b3 = b2;

    // Launch kernels with the bars to modify the foo.
    foo_kernel_launch( b1 );
    foo_kernel_launch( b2 );
    foo_kernel_launch( b3 );

    // Copy the foo back to the host and check the data.
    Foo result = b1.get_foo_on_host();
    EXPECT_EQ( v1 + 3.0, result.get_dbl_on_host() );
    EXPECT_EQ( v2 + 3, result.get_int_on_host() );
}

//---------------------------------------------------------------------------//
// Test argument constructor.
TEST_F(SharedDevicePtrTest, argument_test)
{
    // Test data.
    double v1 = 3.432l;
    int v2 = 4;

    // Make some bars.
    Bar b1( v1, v2 );
    Bar b2 = b1;
    Bar b3 = b2;

    // Launch kernels with the bars to modify the foo.
    foo_kernel_launch( b1 );
    foo_kernel_launch( b2 );
    foo_kernel_launch( b3 );

    // Copy the foo back to the host and check the data.
    Foo result = b1.get_foo_on_host();
    EXPECT_EQ( v1 + 3.0, result.get_dbl_on_host() );
    EXPECT_EQ( v2 + 3, result.get_int_on_host() );
}

//---------------------------------------------------------------------------//
//                        end of tstShared_Device_Ptr.cc
//---------------------------------------------------------------------------//
