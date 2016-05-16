//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Atomic_Add.cuh
 * \author Seth R Johnson
 * \date   Thu Aug 08 08:24:17 2013
 * \brief  Atomic_Add kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Atomic_Add_cuh
#define cuda_utils_Atomic_Add_cuh

#include "Definitions.hh"

namespace cuda
{
//===========================================================================//
/*!
 * \brief Do an atomic add-and-assign
 *
 * \return The old value
 */
template<typename Arch_Switch, typename T>
struct Atomic_Add
{
    typedef Arch_Switch Arch_t;
    typedef T           float_type;

    __device__ float_type operator()(
            float_type* __restrict__  address,
            float_type val);
};

//---------------------------------------------------------------------------//
// HOST EMULATION SPECIALIZATIONS
//---------------------------------------------------------------------------//
//! Specialization on host
template<typename T>
struct Atomic_Add<arch::Host, T>
{
    typedef arch::Host Arch_t;
    typedef T          float_type;

    float_type operator()(float_type* address, float_type val)
    {
        float_type old = *address;
        *address += val;
        return old;
    }
};

#ifdef __NVCC__
//---------------------------------------------------------------------------//
// DEVICE SPECIALIZATIONS
//---------------------------------------------------------------------------//
//! Specialization on device single-precision
template<>
__device__ float Atomic_Add<arch::Device, float>::operator()(
        float* __restrict__  address,
        float val)
{
    return atomicAdd(address, val);
}

//---------------------------------------------------------------------------//
//! Specialization on device double-precision
template<>
__device__ double Atomic_Add<arch::Device, double>::operator()(
        double* __restrict__  address,
        double val)
{
    typedef unsigned long long int ulli;
    ulli* __restrict__ address_as_ull = reinterpret_cast<ulli*>(address);
    ulli old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif // __NVCC__

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_Atomic_Add_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Atomic_Add.cuh
//---------------------------------------------------------------------------//
