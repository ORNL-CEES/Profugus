//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/Exponential.cuh
 * \author Seth R Johnson
 * \date   Fri Aug 16 10:07:03 2013
 * \brief  Exponential kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Exponential_cuh
#define cuda_utils_Exponential_cuh

#ifndef __NVCC__
#include <memory>
#include "Utils/harness/DBC.hh"
#include "Utils/utils/Transport_Exp.hh"
#endif

#include "Definitions.hh"

namespace cuda
{
//===========================================================================//
/*!
 * \brief Do a fast exponential suitable for MOC
 *
 * On the GPU, this uses the \c __expf function, and on the CPU, it uses an
 * exponential table.
 */
template<typename Arch_Switch, typename T>
struct Exponential
{
    typedef Arch_Switch Arch_t;
    typedef T           float_type;

    // For CPU version, set range and precision of exponential
    void initialize(float_type, int) {/* * */}

#ifdef __NVCC__
    __device__ float_type operator()(float_type x);
#endif
};

#ifndef __NVCC__
//---------------------------------------------------------------------------//
// HOST EMULATION SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief The host creates and holds onto an exponential table.
 *
 * It gets passed around as a smart pointer so that kernel calls don't
 * reallocate and recalculate the table (expensive)
 */
template<typename T>
struct Exponential<arch::Host, T>
{
    typedef arch::Host                         Arch_t;
    typedef T                                  float_type;
    typedef profugus::Transport_Exp<float_type> Transport_Exp_t;
    typedef std::shared_ptr<Transport_Exp_t>   SP_Transport_Exp;

    //! Set the range and precision of the exponential
    void initialize(float_type x_end, int intervals)
    {
        Require(!d_table);

        d_table = std::make_shared<Transport_Exp_t>();
        d_table->set_table(x_end, intervals);

        Ensure(d_table);
    }

    //! Calculate the exponential
    float_type operator()(float_type x)
    {
        //Require(d_table);
        return d_table->mexp(-x);
    }

  private:
    SP_Transport_Exp d_table;
};

#else // __NVCC__
//---------------------------------------------------------------------------//
// DEVICE SPECIALIZATIONS
//---------------------------------------------------------------------------//
//! Specialization on device single-precision
template<>
__device__ float Exponential<arch::Device, float>::operator()(float x)
{
    return __expf(x);
}

//---------------------------------------------------------------------------//
//! Specialization on device double-precision
template<>
__device__ double Exponential<arch::Device, double>::operator()(double x)
{
    // Convert from double to float, evaluate, convert from float to double
    return __expf(x);
}

//---------------------------------------------------------------------------//
#endif // __NVCC__

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_Exponential_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Exponential.cuh
//---------------------------------------------------------------------------//
