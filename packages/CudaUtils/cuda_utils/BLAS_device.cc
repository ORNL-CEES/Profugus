//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/BLAS_device.cc
 * \author Seth R Johnson
 * \date   Mon Jul 01 10:42:26 2013
 * \brief  BLAS member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "BLAS.hh"

#include "harness/DBC.hh"

#include "BLAS_device_common.hh"
#include "CudaDBC.hh"

typedef cuda::arch::Device Arch_t;

//---------------------------------------------------------------------------//
// ANONYMOUS NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//
namespace {
/*!
 * \brief Internal function to convert to CUDA enum
 */
cublasOperation_t transpose_convert(const cuda::Transpose t)
{
    cublasOperation_t result = CUBLAS_OP_N;
    switch (t)
    {
        case cuda::NO_TRANS:   result = CUBLAS_OP_N; break;
        case cuda::TRANS:      result = CUBLAS_OP_T; break;
        case cuda::CONJ_TRANS: result = CUBLAS_OP_C; break;
        default:
            CHECK(0);
    }
    return result;
}
} // end anonymous namespace

namespace cuda
{
//---------------------------------------------------------------------------//
// TEMPLATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
template<typename float_type>
void BLAS<Arch_t,float_type>::set_stream(Stream_t& stream)
{
    CudaBlasCall(cublasSetStream(d_blas_handle.handle(), stream.handle()));
}

//---------------------------------------------------------------------------//
// BLAS float SPECIALIZATIONS
//---------------------------------------------------------------------------//
template<>
void BLAS<Arch_t,float>::GEMM(
            Transpose transa, Transpose transb,
            const ord_type m, const ord_type n, const ord_type k,
            const float alpha, const float* a, ord_type lda,
                               const float* b, ord_type ldb,
            const float beta,        float* c, ord_type ldc
            )
{
    cublasOperation_t cutransa = transpose_convert(transa);
    cublasOperation_t cutransb = transpose_convert(transb);

    CudaBlasCall(cublasSgemm(d_blas_handle.handle(), cutransa, cutransb,
            m, n, k,
            &alpha, a, lda,
                    b, ldb,
            &beta,  c, ldc));
}

//---------------------------------------------------------------------------//
template<>
void BLAS<Arch_t,float>::GEMV(
            Transpose transa,
            const ord_type m, const ord_type n,
            const float alpha, const float* a, ord_type lda,
                               const float* x, ord_type incx,
            const float beta,        float* y, ord_type incy
            )
{
    cublasOperation_t cutransa = transpose_convert(transa);

    CudaBlasCall(cublasSgemv(d_blas_handle.handle(), cutransa,
            m, n,
            &alpha, a, lda,
                    x, incx,
            &beta,  y, incy));
}

//---------------------------------------------------------------------------//
// BLAS double SPECIALIZATIONS
//---------------------------------------------------------------------------//
template<>
void BLAS<Arch_t,double>::GEMM(
            Transpose transa, Transpose transb,
            const ord_type m, const ord_type n, const ord_type k,
            const double alpha, const double* a, ord_type lda,
                                const double* b, ord_type ldb,
            const double beta,        double* c, ord_type ldc
            )
{
    cublasOperation_t cutransa = transpose_convert(transa);
    cublasOperation_t cutransb = transpose_convert(transb);

    CudaBlasCall(cublasDgemm(d_blas_handle.handle(), cutransa, cutransb,
            m, n, k,
            &alpha, a, lda,
                    b, ldb,
            &beta,  c, ldc));
}

//---------------------------------------------------------------------------//
template<>
void BLAS<Arch_t,double>::GEMV(
            Transpose transa,
            const ord_type m, const ord_type n,
            const double alpha, const double* a, ord_type lda,
                                const double* x, ord_type incx,
            const double beta,        double* y, ord_type incy
            )
{
    cublasOperation_t cutransa = transpose_convert(transa);

    CudaBlasCall(cublasDgemv(d_blas_handle.handle(), cutransa,
            m, n,
            &alpha, a, lda,
                    x, incx,
            &beta,  y, incy));
}

//---------------------------------------------------------------------------//
// Explicit instantiation
//---------------------------------------------------------------------------//
template class BLAS<Arch_t,float> ;
template class BLAS<Arch_t,double>;
//---------------------------------------------------------------------------//

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of BLAS_device.cc
//---------------------------------------------------------------------------//
