//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstTexture_Vector.cc
 * \author Seth R Johnson
 * \date   Fri Sep 20 15:32:41 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Texture_Vector.hh"

#include <vector>
#include "Utils/gtest/utils_gtest.hh"
#include "Utils/utils/View_Field.hh"

#include "../Hardware.hh"
#include "../Host_Vector.hh"
#include "../Device_Vector.hh"

#include "Texture_Vector_Test_Kernel.cuh"

using cuda::Texture_Vector;

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template <typename Arch_Switch, typename T>
struct Test_Traits
{
    typedef Arch_Switch Arch_t;
    typedef T           value_type;
};

template<typename Test_T>
class TextureVectorTest : public ::testing::Test
{
  protected:
    typedef typename Test_T::Arch_t     Arch_t;
    typedef typename Test_T::value_type value_type;

    // Get device if applicable
    void SetUp()
    {
        typedef cuda::Hardware<Arch_t> Hardware_t;
        if (Hardware_t::valid_device_exists() && !Hardware_t::have_acquired())
        {
            cout << "Acquiring device..." << endl;
            Hardware_t::acquire();
        }
    }
};

typedef Test_Traits<cuda::arch::Host, int>    TT_HI;
typedef Test_Traits<cuda::arch::Host, float>  TT_HF;
typedef Test_Traits<cuda::arch::Host, double> TT_HD;
#ifdef USE_CUDA
typedef Test_Traits<cuda::arch::Device, int>    TT_DI;
typedef Test_Traits<cuda::arch::Device, float>  TT_DF;
typedef Test_Traits<cuda::arch::Device, double> TT_DD;
// instantiate both host and device code
typedef ::testing::Types<TT_HI, TT_HF, TT_HD, TT_DI, TT_DF, TT_DD> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<TT_HI, TT_HF, TT_HD> ArchTypes;
#endif

TYPED_TEST_CASE(TextureVectorTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TYPED_TEST(TextureVectorTest, execute)
{
    typedef typename TypeParam::Arch_t     Arch_t;
    typedef typename TypeParam::value_type value_type;

    typedef cuda::Hardware<Arch_t>                   Hardware_t;
    typedef cuda::Device_Vector<Arch_t, value_type>  Device_Vector_t;
    typedef cuda::Texture_Vector<Arch_t, value_type> Texture_Vector_t;

    if (!Hardware_t::valid_device_exists())
        SKIP_TEST("No valid device exists.");

    // >>> CREATE HOST VECTOR
    std::vector<value_type> original(1024);

    for (std::size_t i = 0; i < original.size(); ++i)
        original[i] = 2 * i + 1;

    // >>> CREATE TEXTURE VECTOR
    Texture_Vector_t texture_vec(original.size());
    EXPECT_EQ(original.size(), texture_vec.size());
    EXPECT_FALSE(texture_vec.is_initialized());

    texture_vec.assign(profugus::make_view(original));

    // >>> RUN (copy from texture to device)
    Device_Vector_t result(texture_vec.size());

    texture_vector_test(texture_vec, result);

    // >>> CHECK
    std::vector<value_type> host_result(original.size());
    device_to_host(result, profugus::make_view(host_result));

    EXPECT_VEC_EQ(original, host_result);

    // >>> REASSIGN TEXTURE

    for (std::size_t i = 0; i < original.size(); ++i)
        original[i] = i + 3;

    texture_vec.assign(profugus::make_view(original));

    // >>> RUN AGAIN
    Device_Vector_t result2(original.size());

    texture_vector_test(texture_vec, result2);

    // >>> CHECK
    device_to_host(result2, profugus::make_view(host_result));

    EXPECT_VEC_EQ(original, host_result);
}

//---------------------------------------------------------------------------//
//                 end of tstTexture_Vector.cc
//---------------------------------------------------------------------------//
