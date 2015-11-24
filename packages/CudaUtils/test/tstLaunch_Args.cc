//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstLaunch_Args.cc
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Launch_Args_Kernel.hh"
#include "../cuda_utils/Launch_Args.hh"
#include "../cuda_utils/Device_Vector.hh"
#include "../cuda_utils/Host_Vector.hh"

//---------------------------------------------------------------------------//
// TEST FIXTURE
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
class Launch_Args_Test : public testing::Test
{
  protected:
    typedef Arch_Switch                        Arch_t;
};

#ifdef USE_CUDA
// instantiate both host and device code
typedef ::testing::Types<cuda::arch::Host, cuda::arch::Device> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<cuda::arch::Host> ArchTypes;
#endif

TYPED_TEST_CASE(Launch_Args_Test, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(Launch_Args_Test, functor)
{
    typedef TypeParam Arch_t;
    
    // Make launch args.
    cuda::Launch_Args<Arch_t> launch_args;
    launch_args.grid_size = 4;
    launch_args.block_size = 256;

    // Create a functor.
    double value = 2.3;
    int size = launch_args.block_size * launch_args.grid_size;
    Functor<Arch_t> functor( size, value );

    // Call the kernel launch.
    cuda::ParallelLaunch<Arch_t>::launch( functor, launch_args );

    // Get the data from the functor.
    const typename Functor<Arch_t>::Host_Vector_t& host_data =
	functor.get_data();

    // Test the data.
    for ( int i = 0; i < size; ++i )
    {
	EXPECT_EQ( host_data[i], value + i );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstLaunch_Args.cc
//---------------------------------------------------------------------------//
