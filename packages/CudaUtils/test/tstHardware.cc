//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstHardware.cc
 * \author Seth R Johnson
 * \date   Wed Jul 10 09:04:11 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "config.h"
#include "../cuda_utils/Hardware.hh"

//---------------------------------------------------------------------------//
#ifdef USE_CUDA
// Make sure we can acquire the device for real
TEST(Initialize, acquire)
{
    typedef cuda::arch::Device     Arch_t;
    typedef cuda::Hardware<Arch_t> Hardware_t;

    if (!Hardware_t::valid_device_exists())
    {
        SKIP_TEST("CUDA is not installed but no valid GPU exists.");
    }

    ASSERT_FALSE(Hardware_t::have_acquired());
    Hardware_t::acquire();
    ASSERT_TRUE(Hardware_t::have_acquired());
}
#endif

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Arch_T>
class HardwareTest : public ::testing::Test
{
  protected:
    typedef Arch_T Arch_t;
};

// instantiate both host and device code
typedef ::testing::Types<cuda::arch::Device, cuda::arch::Host> ArchTypes;

TYPED_TEST_CASE(HardwareTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(HardwareTest, everything)
{
    typedef typename TestFixture::Arch_t Arch_t;
    typedef cuda::Hardware<Arch_t>       Hardware_t;

    // Note: even if CUDA is disabled or not installed, this test should still
    // compile, and the "valid_device_exists" check should still work.
    if (!Hardware_t::valid_device_exists())
    {
        EXPECT_THROW(Hardware_t::acquire(), ::profugus::assertion);
        SKIP_TEST("No valid device exists, or CUDA is not installed.");
    }

    // Acquire device
    if (!Hardware_t::have_acquired())
    {
        Hardware_t::acquire();
    }
    ASSERT_TRUE(Hardware_t::have_acquired());

    cout << "Number of MPs:     " << Hardware_t::num_multiprocessors() << endl;
    cout << "Cores per mp:      " << Hardware_t::num_cores_per_mp() << endl;
    cout << "Threads per warp:  " << Hardware_t::num_threads_per_warp() << endl;
    cout << "Max threads per MP:" << Hardware_t::max_threads_per_mp() << endl;
    cout << "Const memory:      " << Hardware_t::const_memory() << endl;
    cout << "Max shared memory: " << Hardware_t::shared_memory_per_mp() << endl;
    cout << "Number of cores:   " << Hardware_t::num_total_cores() << endl;
}

//---------------------------------------------------------------------------//
//                        end of tstHardware.cc
//---------------------------------------------------------------------------//
