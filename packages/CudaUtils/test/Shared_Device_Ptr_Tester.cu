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
    *(foo->d_dbl_data) += 1.0;
    *(foo->d_int_data) += 1;
}

//---------------------------------------------------------------------------//
// Kernel launch.
void foo_kernel_launch( Bar& bar )
{
    cuda::Shared_Device_Ptr<Foo> foo = bar.get_foo();
    foo_kernel<<<1,1>>>( foo.get_device_ptr() );
}

//---------------------------------------------------------------------------//
// Foo functions
Foo::Foo( const double d, const int i )
{
    cudaMalloc( (void**) &d_dbl_data, sizeof(double) );
    cudaMalloc( (void**) &d_int_data, sizeof(int) );

    cudaMemcpy( d_dbl_data, &d, sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy( d_int_data, &i, sizeof(int), cudaMemcpyHostToDevice );
}

Foo::~Foo()
{
    cudaFree( d_dbl_data );
    cudaFree( d_int_data );
}

double Foo::get_dbl_on_host() const
{
    double host_value = 0.0;
    cudaMemcpy( 
	&host_value, d_dbl_data, sizeof(double), cudaMemcpyDeviceToHost );
    return host_value;
}

int Foo::get_int_on_host() const
{
    int host_value = 0.0;
    cudaMemcpy( &host_value, d_int_data, sizeof(int), cudaMemcpyDeviceToHost );
    return host_value;
}

//---------------------------------------------------------------------------//
// Bar functions.
Bar::Bar( const double d, const int i )
    : d_foo( cuda::shared_device_ptr<Foo>(d, i) )
{ /* ... */ }

Bar::Bar( const std::shared_ptr<Foo>& foo )
    : d_foo( foo )
{ /* ... */ }

const std::shared_ptr<Foo>& Bar::get_foo_on_host() const
{
    return d_foo.get_host_ptr();
}

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Shared_Device_Ptr_Tester.cu
//---------------------------------------------------------------------------//
