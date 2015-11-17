//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstMulti_Vector.cc
 * \author Seth R Johnson
 * \date   Fri Aug 02 12:18:30 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Multi_Vector.hh"

#include "gtest/utils_gtest.hh"
#include "utils/View_Field.hh"

#include "../Hardware.hh"
#include "../Vector_Traits.hh"
#include "../Device_Vector.hh"

using cuda::Multi_Vector;

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Arch_Switch>
class MultiDeviceVectorTest : public ::testing::Test
{
  protected:
    typedef Arch_Switch                         Arch_t;
    typedef cuda::Hardware<Arch_t>              Hardware_t;
    typedef cuda::Vector_Traits<Arch_t, double> V_Traits_t;

    typedef typename V_Traits_t::Multi_Vector_Float  MDV_t;
    typedef typename V_Traits_t::Device_Vector_Float Device_Vector_t;

    typedef std::vector<double> Vector_t;
};


#ifdef USE_CUDA
// instantiate both host and device code
typedef ::testing::Types<cuda::arch::Device, cuda::arch::Host> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<cuda::arch::Host> ArchTypes;
#endif

TYPED_TEST_CASE(MultiDeviceVectorTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TYPED_TEST(MultiDeviceVectorTest, simple)
{
    typedef typename TestFixture::Hardware_t Hardware_t;

    typedef typename TestFixture::Device_Vector_t Device_Vector_t;
    typedef typename TestFixture::Vector_t        Vector_t;

    // >>> GPU INITIALIZATION

    if (!Hardware_t::valid_device_exists())
        SKIP_TEST("No valid device exists.");

    if (!Hardware_t::have_acquired())
    {
        cout << "Acquiring device..." << endl;
        Hardware_t::acquire();
    }

    // Create data on the host
    typename TestFixture::MDV_t mv(10);

    EXPECT_EQ(10, mv.size());
    for (std::size_t i = 0; i < 10; ++i)
        EXPECT_FALSE(mv.is_initialized(i));

    Vector_t host_data(10, 3.);
    mv.initialize(3, profugus::make_view(host_data));

    EXPECT_TRUE(mv.is_initialized(3));
    const Device_Vector_t& dv = mv[3];
    EXPECT_EQ(10, dv.size());

    host_data.assign(5, 2.);
#ifdef REMEMBER_ON
    EXPECT_THROW(mv.initialize(3, profugus::make_view(host_data)),
            profugus::assertion);
#endif

    mv.initialize(5, profugus::make_view(host_data));
    EXPECT_TRUE(mv.is_initialized(5));

    host_data.assign(5, 1.);
    mv.initialize(6, profugus::make_view(host_data));
    swap(mv[5], mv[6]);

    Vector_t result(5);
    cuda::device_to_host(mv[5], profugus::make_view(result));
    EXPECT_VEC_SOFT_EQ(host_data, result);
}

//---------------------------------------------------------------------------//
//                        end of tstMulti_Vector.cc
//---------------------------------------------------------------------------//
