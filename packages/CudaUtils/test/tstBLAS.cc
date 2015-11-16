//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstBLAS.cc
 * \author Seth R Johnson
 * \date   Thu Jul 11 10:12:20 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../BLAS.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "Utils/utils/View_Field.hh"

#include "../Hardware.hh"
#include "../Device_Vector.hh"

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template <typename Arch_Switch, typename T>
struct Test_Traits
{
    typedef Arch_Switch Arch_t;
    typedef T           float_type;
};

//---------------------------------------------------------------------------//
template<typename Traits_T>
class BLASTest : public ::testing::Test
{

  protected:
    typedef typename Traits_T::Arch_t     Arch_t;
    typedef typename Traits_T::float_type float_type;

    typedef cuda::Hardware<Arch_t>                  Hardware_t;
    typedef cuda::BLAS<Arch_t, float_type>          BLAS_t;
    typedef profugus::const_View_Field<float_type>   const_View_Field_t;
    typedef profugus::View_Field<float_type>         View_Field_t;
    typedef cuda::Device_Vector<Arch_t, float_type> Device_Vector_t;

  protected:
    virtual void SetUp()
    {
        // Initialize device
        if (!Hardware_t::have_acquired())
        {
            std::cout << "Acquiring device..." << std::endl;
            Hardware_t::acquire();
        }
        Insist(Hardware_t::have_acquired(), "Device could not be acquired.");
    }

};
//---------------------------------------------------------------------------//

typedef Test_Traits<cuda::arch::Host, float>  TT_HF;
typedef Test_Traits<cuda::arch::Host, double> TT_HD;
#ifdef USE_CUDA
typedef Test_Traits<cuda::arch::Device, float>  TT_DF;
typedef Test_Traits<cuda::arch::Device, double> TT_DD;
// instantiate both host and device code
typedef ::testing::Types<TT_HF, TT_HD, TT_DF, TT_DD> ArchTypes;
#else
// instantiate host-only code
typedef ::testing::Types<TT_HF, TT_HD> ArchTypes;
#endif

TYPED_TEST_CASE(BLASTest, ArchTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TYPED_TEST(BLASTest, GEMV)
{
    typedef typename TestFixture::float_type         float_type;
    typedef typename TestFixture::Hardware_t         Hardware_t;
    typedef typename TestFixture::BLAS_t             BLAS_t;
    typedef typename TestFixture::const_View_Field_t const_View_Field_t;
    typedef typename TestFixture::Device_Vector_t    Device_Vector_t;

    Device_Vector_t matrix_a_gpu(4 * 4);
    Device_Vector_t vector_x_gpu(    4);
    Device_Vector_t vector_y_gpu(    4);

    // Initialize (row-major)
    float_type mat[4][4] = {
        {0.2,     0.4,     0.6,     0.8,},
        {1.2,     1.4,     1.6,     1.8,},
        {2.2,     2.4,     2.6,     2.8,},
        {3.2,     3.4,     3.6,     3.8,},};
    float_type vec[]   = {1.  ,  1.25,  1.5 ,  1.75};

    // Copy to GPU
    matrix_a_gpu.assign(const_View_Field_t(
                &mat[0][0], &mat[0][0] + 4 * 4));
    vector_x_gpu.assign(const_View_Field_t(
                &vec[0],    &vec[0] + 4));

    // Instantiate BLAS
    BLAS_t blas;
    // do y = 1.25 * a * x + 0. * y
    blas.GEMV(cuda::TRANS, 4, 4,
            1.25, matrix_a_gpu.cdata(), 4,
                  vector_x_gpu.cdata(), 1,
            0.,   vector_y_gpu.data(), 1);

    // Wait until kernel is finished
    Hardware_t::synchronize();

    // Copy back to CPU
    std::vector<float_type> result(4);
    cuda::device_to_host(vector_y_gpu, profugus::make_view(result));

    // Test
    std::vector<double> result_dbl(result.begin(), result.end());
    double expected[] = {3.75 ,  10.625,  17.5  ,  24.375};

    EXPECT_VEC_SOFTEQ(expected, result_dbl, 1.e-5);
}

//---------------------------------------------------------------------------//
TYPED_TEST(BLASTest, GEMM)
{
    typedef typename TestFixture::float_type         float_type;
    typedef typename TestFixture::Hardware_t         Hardware_t;
    typedef typename TestFixture::BLAS_t             BLAS_t;
    typedef typename TestFixture::const_View_Field_t const_View_Field_t;
    typedef typename TestFixture::Device_Vector_t    Device_Vector_t;

    Device_Vector_t matrix_a_gpu(4 * 4);
    Device_Vector_t matrix_b_gpu(4 * 4);
    Device_Vector_t matrix_c_gpu(4 * 4);

    // Initialize (row-major)
    float_type mata[4][4] = {
        {0.2,     0.4,     0.6,     0.8,},
        {1.2,     1.4,     1.6,     1.8,},
        {2.2,     2.4,     2.6,     2.8,},
        {3.2,     3.4,     3.6,     3.8,},};
    float_type matb[4][4] = {
        {1.  ,  1.5 ,  2.  ,  2.5 },
        {2.  ,  2.75,  3.5 ,  4.25},
        {3.  ,  4.  ,  5.  ,  6.  },
        {4.  ,  5.25,  6.5 ,  7.75},};

    // Copy to GPU
    matrix_a_gpu.assign(const_View_Field_t(
                &mata[0][0], &mata[0][0] + 4 * 4));
    matrix_b_gpu.assign(const_View_Field_t(
                &matb[0][0], &matb[0][0] + 4 * 4));

    // Instantiate BLAS
    BLAS_t blas;
    // do y = 1.25 * a * x + 0. * y
    blas.GEMM(cuda::TRANS, cuda::TRANS,
            4, 4, 4,
            1.25, matrix_a_gpu.cdata(), 4,
                  matrix_b_gpu.cdata(), 4,
            0.,   matrix_c_gpu.data(), 4);

    // Wait until kernel is finished
    Hardware_t::synchronize();

    // Copy back to CPU
    std::vector<float_type> result(4 * 4);
    cuda::device_to_host(matrix_c_gpu, profugus::make_view(result));

    // Test
    std::vector<double> result_dbl(result.begin(), result.end());
    // Linearized column-major result
    double expected[] = {
        7.5   ,  20.   ,  32.5  ,  45.   ,  10.   ,  26.875,  43.75 ,
        60.625,  12.5  ,  33.75 ,  55.   ,  76.25 ,  15.   ,  40.625,
        66.25 ,  91.875};

    EXPECT_VEC_SOFTEQ(expected, result_dbl, 1.e-5);
}

//---------------------------------------------------------------------------//
//                        end of tstBLAS.cc
//---------------------------------------------------------------------------//
