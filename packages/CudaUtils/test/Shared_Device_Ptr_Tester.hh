//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Shared_Device_Ptr_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Shared_Device_Ptr_Tester_hh
#define cuda_utils_test_Shared_Device_Ptr_Tester_hh

#include <memory>

#include "../cuda_utils/Shared_Device_Ptr.hh"

//---------------------------------------------------------------------------//
// Helper classes to test deep/shallow copy.
//---------------------------------------------------------------------------//
class Foo
{
  public:

    Foo() { /* ... */ }

    Foo( const double d, const int i );

    ~Foo();

    Foo( const Foo& foo ) = default;
    Foo& operator=( const Foo& foo ) = default;

    double get_dbl_on_host() const;

    int get_int_on_host() const;

  public:

    // DEVICE DATA
    double* d_dbl_data;
    int* d_int_data;
};

//---------------------------------------------------------------------------//
class Bar
{
  public:
    Bar( const double d, const int i );

    Bar( const std::shared_ptr<Foo>& foo );

    Foo get_foo_on_host() const;

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
