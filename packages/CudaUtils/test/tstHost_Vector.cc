//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstHost_Vector.cc
 * \author Seth R Johnson
 * \date   Mon Aug 12 10:29:26 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Host_Vector.hh"

#include "gtest/utils_gtest.hh"

#include "../Device_Vector.hh"
#include "../Hardware.hh"
#include "utils/View_Field.hh"

#include <../config.h>
#include "Polyglot_Kernel.cuh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class HostVectorTest : public ::testing::Test
{
  protected:
    typedef cuda::arch::Device                Arch_t;
    typedef cuda::Device_Vector<Arch_t,float> Device_Vector_t;

    typedef std::vector<float>               Vector_t;
    typedef cuda::Host_Vector<float>         Host_Vector_t;
    typedef profugus::const_View_Field<float> const_View_Field_t;
    typedef profugus::View_Field<float>       View_Field_t;

  protected:
    void SetUp()
    {
#ifdef USE_CUDA
        typedef cuda::Hardware<Arch_t> Hardware_t;
        // Initialize device
        if (!Hardware_t::have_acquired())
        {
            std::cout << "Acquiring device..." << std::flush;
            Hardware_t::acquire();
            std::cout << "done." << std::endl;
        }
        Insist(Hardware_t::have_acquired(), "Device could not be acquired.");
#endif

        // Add values to the vector
        for (std::size_t i = 0; i < 63; ++i)
        {
            original.push_back(i * i);
        }
    }

  protected:
    Vector_t original;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HostVectorTest, accessors)
{
    Host_Vector_t hv(original.size(), 1.23f);

    ASSERT_EQ(63, hv.size());
    ASSERT_TRUE(!hv.empty());
    EXPECT_FALSE(hv.is_mapped());
    EXPECT_FALSE(hv.is_write_combined());

    // Check default value and iterators
    for (Host_Vector_t::const_iterator it = hv.begin(),
            end_it = hv.end();
            it != end_it;
            ++it)
    {
        EXPECT_FLOAT_EQ(1.23f, *it);
    }

    // Copy data
    hv.assign(profugus::make_view(original));

    // Check equivalence and bracket accessors
    for (std::size_t i = 0; i < 63; ++i)
    {
        EXPECT_FLOAT_EQ(original[i], hv[i]);
    }
}

//---------------------------------------------------------------------------//
// Mapped memory is only supported when CUDA is enabled
TEST_F(HostVectorTest, mapped_memory)
{
#ifdef USE_CUDA
    // Create a "mapped memory, write-only" vector where assigning to the host
    // memory will automagically be copied to the GPU
    Host_Vector_t hv(original.size(), 0.f,
            cuda::alloc::MAPPED_WRITE_COMBINED);

    // Create blank destination vector
    Device_Vector_t gpu_out(original.size());

    // Copy data
    hv.assign(profugus::make_view(original));

    // Call kernel; the host wrapper uses the .data() accessor to extract the
    // GPU pointer
    polyglot_copy(hv, gpu_out);

    // Copy from device vector to new vector
    Vector_t result(original.size());
    device_to_host(gpu_out, profugus::make_view(result));

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], result[i]) << "Failure at index " << i;
    }
#else
    // No cuda, no mapped memory; Insist should be raised
    EXPECT_THROW({Host_Vector_t hv(original.size(), 0.f,
            cuda::alloc::MAPPED_WRITE_COMBINED);}, profugus::assertion);
#endif
}
//---------------------------------------------------------------------------//
//                        end of tstHost_Vector.cc
//---------------------------------------------------------------------------//
