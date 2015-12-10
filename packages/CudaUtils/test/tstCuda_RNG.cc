//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstCuda_RNG.cc
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Run_Cuda_RNG.hh"
#include "../cuda_utils/Cuda_RNG.hh"
#include "../cuda_utils/Device_Vector.hh"
#include "../cuda_utils/Host_Vector.hh"

//---------------------------------------------------------------------------//
// TEST FIXTURE
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
class Cuda_RNG_Test : public testing::Test
{
};

#ifdef USE_CUDA
// instantiate both host and device code
typedef ::testing::Types<cuda::arch::Host, cuda::arch::Device> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<cuda::arch::Host> ArchTypes;
#endif

TYPED_TEST_CASE(Cuda_RNG_Test, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(Cuda_RNG_Test, functor)
{
    typedef TypeParam Arch_t;

    // Size of test.
    int num_ran = 100;

    // Create a cpu rng to get started.
    cuda::Cuda_Global_RNG_Control::initialize( 393243 );
    profugus::RNG cpu_rng = cuda::Cuda_Global_RNG_Control::d_rng_control.rng(0);
    
    // Create seeds to initialize the vector with.
    cuda::Host_Vector<int> host_seeds( num_ran );
    for ( int i = 0; i < num_ran; ++i )
    {
	host_seeds[i] = cpu_rng.ran();
    }

    // Create a vector or rngs.
    cuda::Host_Vector<cuda::Cuda_RNG> host_rngs( num_ran );

    // Create a vector to fill with random numbers.
    cuda::Host_Vector<double> host_data( num_ran );

    // Fill the vector.
    run_cuda_rng<Arch_t>( host_seeds, host_rng, host_data );

    // Test the data.
    for ( int i = 0; i < host_data.size(); ++i )
    {
	std::cout << host_data[i] << std::endl;
    }
}

//---------------------------------------------------------------------------//
//                        end of tstCuda_RNG.cc
//---------------------------------------------------------------------------//
