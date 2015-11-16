//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/BLAS_host.cc
 * \author Seth R Johnson
 * \date   Wed Dec 11 20:33:47 2013
 * \brief  BLAS member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "BLAS.hh"

#include "Teuchos_BLAS.hpp"

#include "Utils/harness/DBC.hh"

typedef cuda::arch::Host Arch_t;

//---------------------------------------------------------------------------//
// ANONYMOUS NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//
namespace {
/*!
 * \brief Internal function to convert to CUDA enum
 */
inline Teuchos::ETransp transpose_convert(const cuda::Transpose t)
{
    Teuchos::ETransp result = Teuchos::NO_TRANS;
    switch (t)
    {
        case cuda::NO_TRANS:   result = Teuchos::NO_TRANS;   break;
        case cuda::TRANS:      result = Teuchos::TRANS;      break;
        case cuda::CONJ_TRANS: result = Teuchos::CONJ_TRANS; break;
        default:
            Check(0);
    }
    return result;
}
} // end anonymous namespace

namespace cuda
{

//---------------------------------------------------------------------------//
// BLAS
//---------------------------------------------------------------------------//
template<typename T>
void BLAS<Arch_t,T>::GEMM(
            Transpose transa, Transpose transb,
            const ord_type m, const ord_type n, const ord_type k,
            const float_type alpha, const float_type* a, ord_type lda,
                                    const float_type* b, ord_type ldb,
            const float_type beta,        float_type* c, ord_type ldc
            )
{
    Teuchos::ETransp ttransa = transpose_convert(transa);
    Teuchos::ETransp ttransb = transpose_convert(transb);

    Teuchos_Blas_t blas;

    blas.GEMM(ttransa, ttransb,
            m, n, k,
            alpha, a, lda,
                   b, ldb,
            beta,  c, ldc);
}

//---------------------------------------------------------------------------//
template<typename T>
void BLAS<Arch_t,T>::GEMV(
            Transpose transa,
            const ord_type m, const ord_type n,
            const float_type alpha, const float_type* a, ord_type lda,
                                    const float_type* x, ord_type incx,
            const float_type beta,        float_type* y, ord_type incy
            )
{
    Teuchos::ETransp ttransa = transpose_convert(transa);

    Teuchos_Blas_t blas;

    blas.GEMV(ttransa,
            m, n,
            alpha, a, lda,
                   x, incx,
            beta,  y, incy);
}

//---------------------------------------------------------------------------//
// Explicit instantiation
//---------------------------------------------------------------------------//
template class BLAS<Arch_t,float> ;
template class BLAS<Arch_t,double>;
//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of BLAS_host.cc
//---------------------------------------------------------------------------//
