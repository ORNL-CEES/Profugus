//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstAtomic_Add.cc
 * \author Seth R Johnson
 * \date   Thu Aug 15 11:22:01 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Utils/gtest/utils_gtest.hh"

#include "Atomic_Add_Kernel_Data.hh"
#include "Atomic_Add_Kernel.cuh"
#include "../Hardware.hh"
#include "../Host_Vector.hh"

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Vector_Traits>
class LockTest : public ::testing::Test
{
};

typedef cuda::Vector_Traits<cuda::arch::Host, float>  VT_HF;
typedef cuda::Vector_Traits<cuda::arch::Host, double> VT_HD;
#ifdef USE_CUDA
typedef cuda::Vector_Traits<cuda::arch::Device, float>  VT_DF;
typedef cuda::Vector_Traits<cuda::arch::Device, double> VT_DD;
// instantiate both host and device code
typedef ::testing::Types<VT_HF, VT_HD, VT_DF, VT_DD> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<VT_HF> ArchTypes;
#endif

TYPED_TEST_CASE(LockTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(LockTest, execute)
{
    typedef TypeParam                            Vector_Traits_t;
    typedef typename Vector_Traits_t::Arch_t     Arch_t;
    typedef typename Vector_Traits_t::float_type float_type;
    typedef cuda::Hardware<Arch_t>               Hardware_t;

    typedef cuda::Atomic_Add_Kernel_Data<Arch_t, float_type> Kernel_Data_t;

    typedef typename Vector_Traits_t::Host_Vector_Float Host_Vector_Float;

    if (!Hardware_t::valid_device_exists())
        SKIP_TEST("No valid device exists.");

    if (!Hardware_t::have_acquired())
    {
        cout << "Acquiring device..." << endl;
        Hardware_t::acquire();
    }

    // Allocate kernel data
    Kernel_Data_t data;

    // Initialize kernel data
    {
        Host_Vector_Float zeroes(1, static_cast<float_type>(0));
        data.output.assign(zeroes);
    }

    // Set number of threads and blocks to reasonable GPU/CPU numbers
    data.launch_args.block_size = Hardware_t::num_cores_per_mp();
    data.launch_args.grid_size  = Hardware_t::num_multiprocessors();

    // Number of increments to perform
    data.num_increments = 10000;

    // Do the kernel call
    cuda::atomic_add_test(data);

    // Check results
    Host_Vector_Float cpu_result(data.output.size());
    cuda::device_to_host(data.output, cpu_result);

    EXPECT_EQ(1, cpu_result.size());
    EXPECT_EQ(data.num_increments, cpu_result[0]);
}

//---------------------------------------------------------------------------//
//                        end of tstAtomic_Add.cc
//---------------------------------------------------------------------------//
