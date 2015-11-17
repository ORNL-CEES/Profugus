//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstExponential.cc
 * \author Seth R Johnson
 * \date   Fri Aug 16 10:22:07 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Exp_Kernel.cuh"
#include "../Vector_Traits.hh"
#include "../Host_Vector.hh"
#include "../Device_Vector.hh"
#include "../Hardware.hh"
#include "../Exponential.cuh"

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Vector_Traits>
class ExpTest : public ::testing::Test
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

TYPED_TEST_CASE(ExpTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(ExpTest, execute)
{
    typedef TypeParam                            Vector_Traits_t;
    typedef typename Vector_Traits_t::Arch_t     Arch_t;
    typedef typename Vector_Traits_t::float_type float_type;
    typedef cuda::Hardware<Arch_t>               Hardware_t;

    typedef cuda::Exponential<Arch_t, float_type> Exponential_t;

    typedef typename Vector_Traits_t::Device_Vector_Float Device_Vector_Float;
    typedef typename Vector_Traits_t::Host_Vector_Float   Host_Vector_Float;

    if (!Hardware_t::valid_device_exists())
        SKIP_TEST("No valid device exists.");

    if (!Hardware_t::have_acquired())
    {
        cout << "Acquiring device..." << endl;
        Hardware_t::acquire();
    }

    std::vector<double> expected(1024);

    // Create exponential table if host; null-op if GPU
    Exponential_t exp;
    exp.initialize(static_cast<float_type>(expected.size()) / 100, 10000);

    // Allocate in/out/result data
    Device_Vector_Float data(expected.size());

    // Initialize kernel data
    {
        Host_Vector_Float host(data.size());
        for (std::size_t i = 0; i < host.size(); ++i)
        {
            host[i] = -static_cast<float_type>(i) / 100;
            expected[i] = std::exp(host[i]);
        }
        data.assign(host);
    }

    // Do the kernel call
    cuda::exp_test(data, exp);

    // Check results
    Host_Vector_Float cpu_result(data.size());
    cuda::device_to_host(data, cpu_result);

    std::vector<double> actual(cpu_result.begin(), cpu_result.end());
    EXPECT_VEC_SOFTEQ(expected, actual, 1.e-6);
}

//---------------------------------------------------------------------------//
//                        end of tstExponential.cc
//---------------------------------------------------------------------------//
