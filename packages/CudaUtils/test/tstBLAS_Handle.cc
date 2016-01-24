//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstBLAS_Handle.cc
 * \author Seth R Johnson
 * \date   Wed Dec 11 21:01:29 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../cuda_utils/BLAS_Handle.hh"

#include "gtest/utils_gtest.hh"
#include "../cuda_utils/Hardware.hh"

using cuda::BLAS_Handle;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class BLASTest : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
        typedef cuda::arch::Device     Arch_t;
        typedef cuda::Hardware<Arch_t> Hardware_t;

        // Initialize device
        if (!Hardware_t::have_acquired())
        {
            std::cout << "Acquiring device..." << std::endl;
            Hardware_t::acquire();
        }
        INSIST(Hardware_t::have_acquired(), "Device could not be acquired.");
    }
};

//---------------------------------------------------------------------------//
TEST_F(BLASTest, handle_count)
{
    using cuda::BLAS_Handle;

    ASSERT_EQ(0, BLAS_Handle::handle_count());
    {
        BLAS_Handle a;
        EXPECT_EQ(1, BLAS_Handle::handle_count());
        BLAS_Handle b;
        EXPECT_EQ(2, BLAS_Handle::handle_count());
    }
    EXPECT_EQ(0, BLAS_Handle::handle_count());
}

//---------------------------------------------------------------------------//
//                 end of tstBLAS_Handle.cc
//---------------------------------------------------------------------------//
