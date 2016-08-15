//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstDevice_Vector.cc
 * \author Seth R Johnson
 * \date   Thu Aug 01 13:12:15 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../cuda_utils/Device_Vector.hh"

#include <vector>
#include "gtest/utils_gtest.hh"
#include "utils/View_Field.hh"

#include "../cuda_utils/Hardware.hh"
#include "../cuda_utils/Host_Vector.hh"
#include "Polyglot_Kernel.cuh"

#include <type_traits>

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template<typename Arch_Switch>
class DeviceVectorTest : public ::testing::Test
{
  protected:
    typedef Arch_Switch                       Arch_t;
    typedef cuda_utils::Device_Vector<Arch_t,float> Device_Vector_t;

    typedef std::vector<float>               Vector_t;
    typedef profugus::const_View_Field<float> const_View_Field_t;
    typedef profugus::View_Field<float> View_Field_t;
    typedef cuda_utils::Host_Vector<float>         Host_Vector_t;

  protected:
    void SetUp()
    {
        typedef cuda_utils::Hardware<Arch_t> Hardware_t;

        // Initialize device
        if (!Hardware_t::have_acquired())
        {
            std::cout << "Acquiring device..." << std::flush;
            Hardware_t::acquire();
            std::cout << "done." << std::endl;
        }
        INSIST(Hardware_t::have_acquired(), "Device could not be acquired.");

        // Add values to the vector
        for (std::size_t i = 0; i < 63; ++i)
        {
            original.push_back(i * i);
        }
    }

    void takes_a_const(const float*)
    {
        /* * */
    }

    void takes_a_mutable(float*)
    {
        /* * */
    }


  protected:
    Vector_t original;
};

#ifdef USE_CUDA
// instantiate both host and device code
typedef ::testing::Types<cuda_utils::arch::Host, cuda_utils::arch::Device> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<cuda_utils::arch::Host> ArchTypes;
#endif

TYPED_TEST_CASE(DeviceVectorTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, initialize)
{
    typedef typename TestFixture::const_View_Field_t const_View_Field_t;
    typedef typename TestFixture::Device_Vector_t    Device_Vector_t;
    typedef typename TestFixture::Vector_t           Vector_t;

    const Vector_t& original = this->original;

    // Create copy of data on GPU
    const_View_Field_t input_view = profugus::make_view(original);
    Device_Vector_t gpu_in(input_view);
    // Create blank vector
    Device_Vector_t gpu_out(original.size());

    // Call kernel to copy data
    polyglot_copy(gpu_in, gpu_out);

    // Copy from device vector to new vector
    Vector_t result(original.size());
    device_to_host(gpu_out, profugus::make_view(result));

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], result[i]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, assign)
{
    typedef typename TestFixture::Device_Vector_t Device_Vector_t;
    typedef typename TestFixture::Vector_t        Vector_t;

    const Vector_t& original = this->original;

    // Create blank vectors
    Device_Vector_t gpu_in(original.size());

    // Assign
    gpu_in.assign(profugus::make_view(original));

    // Copy from device vector to new vector
    Vector_t result(original.size());
    device_to_host(gpu_in, profugus::make_view(result));

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], result[i]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, assign_device)
{
    typedef typename TestFixture::Device_Vector_t Device_Vector_t;
    typedef typename TestFixture::Vector_t        Vector_t;

    const Vector_t& original = this->original;

    // Create blank vectors
    Device_Vector_t gpu_in(original.size());
    Device_Vector_t gpu_out(original.size());

    // Assign
    gpu_in.assign(profugus::make_view(original));
    gpu_out.assign(gpu_in);


    // Copy from device vector to new vector
    Vector_t result(original.size());
    device_to_host(gpu_out, profugus::make_view(result));

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], result[i]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, write_combined_memory)
{
    typedef typename TestFixture::Arch_t          Arch_t;
    typedef typename TestFixture::Device_Vector_t Device_Vector_t;
    typedef typename TestFixture::Host_Vector_t   Host_Vector_t;
    typedef typename TestFixture::Vector_t        Vector_t;

    using namespace cuda_utils::alloc;

    const Vector_t& original = this->original;

    // Create blank vectors
    Device_Vector_t gpu_in(original.size());

    // If we're using host emulation with CUDA installed, we can't use
    // write-combined memory: hence adjust_alloc_flag<Arch_t>
    // Assign data to write-only host vector
    Host_Vector_t host_data(original.size(), 0.f,
            cuda_utils::adjust_alloc_flag<Arch_t>(WRITE_COMBINED));
    std::copy(original.begin(), original.end(), host_data.begin());

    // Copy to device vector
    gpu_in.assign(host_data);

    // Copy from device vector to new vector
    Vector_t result(original.size());
    device_to_host(gpu_in, profugus::make_view(result));

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], result[i]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, mapped_memory)
{
#ifndef USE_CUDA
    SKIP_TEST("Can't test mapped memory because CUDA is disabled");
#endif
    typedef typename TestFixture::Device_Vector_t Device_Vector_t;
    typedef typename TestFixture::Host_Vector_t   Host_Vector_t;
    typedef typename TestFixture::Vector_t        Vector_t;

    const Vector_t& original = this->original;

    Host_Vector_t a(original.size(), 0, cuda_utils::alloc::MAPPED);
    // Assign to mapped memory
    std::copy(original.begin(), original.end(), a.begin());

    // Create device vector
    Device_Vector_t b(original.size());

    // Copy to device vector (device-to-device hopefully)
    b.assign(a);

    // Create host vector (mapped)
    Host_Vector_t c(original.size(), 0, cuda_utils::alloc::MAPPED);
    device_to_host(b, c);

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], c[i]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, copy_constructor)
{
    typedef typename TestFixture::const_View_Field_t const_View_Field_t;
    typedef typename TestFixture::Device_Vector_t    Device_Vector_t;
    typedef typename TestFixture::Vector_t           Vector_t;

    const Vector_t& original = this->original;

    // Create blank vectors
    const_View_Field_t input_view = profugus::make_view(original);
    Device_Vector_t gpu_in(input_view);

    // Copy constructor
    Device_Vector_t dupe(gpu_in);

    // Copy from device vector to new vector
    Vector_t result(original.size());
    device_to_host(dupe, profugus::make_view(result));

    // Check result
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_FLOAT_EQ(original[i], result[i]);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(DeviceVectorTest, initialized)
{
    typedef typename TestFixture::const_View_Field_t const_View_Field_t;
    typedef typename TestFixture::Device_Vector_t    Device_Vector_t;
    typedef typename TestFixture::Vector_t           Vector_t;

    const Vector_t& original = this->original;

    // Create blank vectors
    const_View_Field_t input_view = profugus::make_view(original);
    Device_Vector_t gpu_a(input_view);
    Device_Vector_t gpu_b(original.size());

    EXPECT_TRUE(gpu_a.is_initialized());
    EXPECT_FALSE(gpu_b.is_initialized());

    // Assign
    gpu_b.assign(profugus::make_view(original));
    EXPECT_TRUE(gpu_b.is_initialized());
}

TYPED_TEST(DeviceVectorTest, accessors)
{
    typedef typename TestFixture::Device_Vector_t Device_Vector_t;
    typedef typename TestFixture::Vector_t        Vector_t;

    const Vector_t& original = this->original;

    Device_Vector_t gpu_b(original.size());

#ifdef REQUIRE_ON
    EXPECT_THROW(this->takes_a_const(gpu_b.cdata()), profugus::assertion);

    const Device_Vector_t& gpu_b_const = gpu_b;
    EXPECT_THROW(gpu_b_const.data(), profugus::assertion);
#endif

    // This should change it to initialized
    EXPECT_NO_THROW(this->takes_a_mutable(gpu_b.data()));
    EXPECT_TRUE(gpu_b.is_initialized());
}

//---------------------------------------------------------------------------//
TYPED_TEST( DeviceVectorTest, device_api )
{
    typedef typename TestFixture::const_View_Field_t const_View_Field_t;
    typedef typename TestFixture::View_Field_t       View_Field_t;
    typedef typename TestFixture::Device_Vector_t    Device_Vector_t;
    typedef typename TestFixture::Vector_t           Vector_t;

    // This test is for device-only vectors.
    if ( std::is_same<typename TestFixture::Arch_t,cuda_utils::arch::Device>::value )
    {
	const Vector_t& original = this->original;

	// Create copy of data on GPU
	const_View_Field_t input_view = profugus::make_view(original);
	Device_Vector_t gpu_in(input_view);

	// Create blank vector
	Vector_t result(original.size());
	View_Field_t output_view = profugus::make_view(result);
	Device_Vector_t gpu_out(output_view);

	// Call kernel to copy data
	polyglot_copy_vector(gpu_in, gpu_out);

	// Copy from device vector to new vector
	device_to_host(gpu_out, profugus::make_view(result));

	// Check result
	for (std::size_t i = 0; i < original.size(); ++i)
	{
	    EXPECT_FLOAT_EQ(original[i], result[i]);
	}
    }
}

//---------------------------------------------------------------------------//
//                        end of tstDevice_Vector.cc
//---------------------------------------------------------------------------//
