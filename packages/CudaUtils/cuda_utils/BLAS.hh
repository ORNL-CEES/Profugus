//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/BLAS.hh
 * \author Seth R Johnson
 * \date   Mon Jul 01 10:42:48 2013
 * \brief  BLAS class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_BLAS_hh
#define cuda_utils_BLAS_hh

#include "config.h"

#include "harness/DBC.hh"

#include "Definitions.hh"
#include "Stream.hh"
#ifdef USE_CUDA
#include "BLAS_Handle.hh"
#endif

// Declare Teuchos classes for Host code
namespace Teuchos
{
template<typename OrdinalType, typename ScalarType> class BLAS;
}

namespace cuda_utils
{
//===========================================================================//
/*!
 * \enum Transpose
 * \brief CUDA interface for BLAS, to be used from host code with device ptrs
 */
enum Transpose
{
    NO_TRANS = 0, //! Not transposed
    TRANS,        //! Transposed
    CONJ_TRANS    //! Conjugate transposed
};

//===========================================================================//
/*!
 * \class BLAS
 * \brief Interface for BLAS compatible with CPUs and GPUs
 */
//===========================================================================//
template <typename Arch_T, typename T>
class BLAS
{
  private:
    // Only use specializations!
    BLAS();
};

//===========================================================================//
// HOST SPECIALIZATION
//===========================================================================//
template<class T>
class BLAS<arch::Host, T>
{
    typedef BLAS<arch::Host, T> This;
  public:
    //! Architecture type
    typedef cuda_utils::arch::Host Arch_t;
    //! Ordinal type
    typedef int ord_type;
    //! Floating point type
    typedef T float_type;
    //! Stream type
    typedef Stream<Arch_t> Stream_t;

  private:
    // Teuchos BLAS interface type.
    typedef Teuchos::BLAS<int, float_type> Teuchos_Blas_t;

  public:
    //! Initialization needs nothing
    BLAS() {/* * */}

    //! Set the stream for the next kernel call (null-op on CPU code)
    void set_stream(Stream_t&) { /* * */ }

    // Matrix-matrix multiply
    void GEMM(
            Transpose transa, Transpose transb,
            const ord_type m, const ord_type n, const ord_type k,
            const float_type alpha, const float_type* a, ord_type lda,
                                    const float_type* b, ord_type ldb,
            const float_type beta,        float_type* c, ord_type ldc
            );

    // Matrix-vector multiply
    void GEMV(
            Transpose transa,
            const ord_type m, const ord_type n,
            const float_type alpha, const float_type* a, ord_type lda,
                                    const float_type* x, ord_type incx,
            const float_type beta,        float_type* y, ord_type incy
            );
};

#ifdef USE_CUDA
//===========================================================================//
// DEVICE SPECIALIZATION
//===========================================================================//
template<class T>
class BLAS<arch::Device, T>
{
    typedef BLAS<arch::Device, T> This;
  public:
    //! Architecture type
    typedef cuda_utils::arch::Device Arch_t;
    //! Ordinal type
    typedef int ord_type;
    //! Floating point type
    typedef T float_type;
    //! Stream type
    typedef Stream<Arch_t> Stream_t;

  private:
    //! Wrapper for creation/destruction of BLAS calls
    BLAS_Handle d_blas_handle;

  public:
    //! Initialize with automatic handle instantiation
    BLAS() { /* * */ }

    //! Set up the stream for the next matrix call
    void set_stream(Stream_t& stream);

    // Matrix-matrix multiply
    void GEMM(
            Transpose transa, Transpose transb,
            const ord_type m, const ord_type n, const ord_type k,
            const float_type alpha, const float_type* a, ord_type lda,
                                    const float_type* b, ord_type ldb,
            const float_type beta,        float_type* c, ord_type ldc
            );

    // Matrix-vector multiply
    void GEMV(
            Transpose transa,
            const ord_type m, const ord_type n,
            const float_type alpha, const float_type* a, ord_type lda,
                                    const float_type* x, ord_type incx,
            const float_type beta,        float_type* y, ord_type incy
            );
};
#endif // USE_CUDA
//===========================================================================//

} // end namespace cuda_utils

#endif // cuda_utils_BLAS_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/BLAS.hh
//---------------------------------------------------------------------------//
