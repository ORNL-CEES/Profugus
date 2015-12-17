//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Shared_Device_Ptr_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Shared_Device_Ptr_Tester_hh
#define cuda_utils_test_Shared_Device_Ptr_Tester_hh

#include "../cuda_utils/Shared_Device_Ptr.hh"

//---------------------------------------------------------------------------//
// Helper classes.
//---------------------------------------------------------------------------//
class Foo
{
  public:

    Foo() { /* ... */ }

    Foo( const double d, const int i )
	: d_dbl_data( d )
	, d_int_data( i )
    { /* ... */ }

  public:
    double d_dbl_data;
    int d_int_data;
};

//---------------------------------------------------------------------------//
class Bar
{
  public:
    Bar( const double d, const int i );

    Bar( const Foo& foo );

    Foo get_foo_on_host();

    cuda::Shared_Device_Ptr<Foo> get_foo() { return d_foo; }

  private:
    cuda::Shared_Device_Ptr<Foo> d_foo;
};

//---------------------------------------------------------------------------//
// Kernel Launch
void foo_kernel_launch( Bar& bar );

//---------------------------------------------------------------------------//
#endif // cuda_utils_test_Shared_Device_Ptr_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Shared_Device_Ptr_Tester.cuh
//---------------------------------------------------------------------------//
