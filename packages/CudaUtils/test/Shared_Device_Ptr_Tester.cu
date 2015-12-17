//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Shared_Device_Ptr_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Shared_Device_Ptr_Tester.hh"

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// Test kernel.
__global__ void foo_kernel( Foo* foo )
{
    foo->d_dbl_data += 1.0;
    foo->d_int_data += 1;
}

//---------------------------------------------------------------------------//
// Kernel launch.
void foo_kernel_launch( Bar& bar )
{
    cuda::Shared_Device_Ptr<Foo> foo = bar.get_foo();
    foo_kernel<<<1,1>>>( foo.get_device_ptr() );
}

//---------------------------------------------------------------------------//
// Bar functions.
Bar::Bar( const double d, const int i )
    : d_foo( d, i )
{ /* ... */ }

Bar::Bar( const Foo& foo )
    : d_foo( foo )
{ /* ... */ }

Foo Bar::get_foo_on_host()
{
    Foo foo;
    cudaMemcpy( &foo, d_foo.get_device_ptr(), sizeof(Foo),
		cudaMemcpyDeviceToHost );
    return foo;
}

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Shared_Device_Ptr_Tester.cu
//---------------------------------------------------------------------------//
