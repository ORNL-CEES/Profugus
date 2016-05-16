//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/test/tstLaunch_Args.cc
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Run_Launch_Args.hh"
#include "../cuda_utils/Launch_Args.hh"
#include "../cuda_utils/Device_Vector.hh"
#include "../cuda_utils/Host_Vector.hh"

//---------------------------------------------------------------------------//
// TEST FIXTURE
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
class Launch_Args_Test : public testing::Test
{
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

    // Call the kernel launch.
    std::vector<double> host_data;
    run_launch_args<Arch_t>( host_data );

    double value = 2.3;

    // Test the data.
    for ( int i = 0; i < host_data.size(); ++i )
    {
        EXPECT_EQ( host_data[i], value + i );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstLaunch_Args.cc
//---------------------------------------------------------------------------//
