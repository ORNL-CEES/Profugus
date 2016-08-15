//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstAtomic_Lock.cc
 * \author Seth R Johnson
 * \date   Thu Aug 15 08:19:34 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../cuda_utils/Hardware.hh"
#include "../cuda_utils/Vector_Traits.hh"
#include "../cuda_utils/Host_Vector.hh"

#include "Atomic_Lock_Test_Kernel.cuh"

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
class LockTest : public ::testing::Test
{
  protected:
    typedef Arch_Switch                 Arch_t;
    typedef cuda_utils::Hardware<Arch_t>      Hardware_t;
    typedef cuda_utils::Vector_Traits<Arch_t> Vector_Traits_t;
};

#ifdef USE_CUDA
// instantiate both host and device code
typedef ::testing::Types<cuda_utils::arch::Device, cuda_utils::arch::Host> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<cuda_utils::arch::Host> ArchTypes;
#endif

TYPED_TEST_CASE(LockTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(LockTest, execute)
{
    typedef typename TestFixture::Arch_t          Arch_t;
    typedef typename TestFixture::Hardware_t      Hardware_t;
    typedef typename TestFixture::Vector_Traits_t Vector_Traits_t;

    typedef cuda_utils::Lock_Kernel_Data<Arch_t> Kernel_Data_t;

    typedef typename Vector_Traits_t::Host_Vector_Float Host_Vector_Float;


    if (!Hardware_t::valid_device_exists())
        SKIP_TEST("No valid device exists.");

    if (!Hardware_t::have_acquired())
    {
        cout << "Acquiring device..." << endl;
        Hardware_t::acquire();
    }

    // Create data on the host
    Host_Vector_Float original(1234);
    for (unsigned int i = 0; i < original.size(); ++i)
        original[i] = 0.;

    // Allocate kernel data
    Kernel_Data_t data(original.size());

    // Initialize data
    data.output.assign(original);

    // Set number of threads and blocks to reasonable GPU/CPU numbers
    data.launch_args.set_block_size( Hardware_t::num_cores_per_mp() );
    data.launch_args.set_grid_size(  Hardware_t::num_multiprocessors() );

    // Do the kernel call
    cuda_utils::lock_test(data);

    // Check results
    Host_Vector_Float cpu_result(data.output.size());
    cuda_utils::device_to_host(data.output, cpu_result);

    EXPECT_EQ(original.size(), cpu_result.size());
    for (unsigned int i = 0; i < original.size(); ++i)
        EXPECT_FLOAT_EQ(data.launch_args.grid_size(), cpu_result[i]);
}

//---------------------------------------------------------------------------//
//                        end of tstAtomic_Lock.cc
//---------------------------------------------------------------------------//
