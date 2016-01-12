//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/Atomic_Add.cuh
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
__inline__ __device__ double Atomic_Add<arch::Device, double>::operator()(
        double* __restrict__  address,
        double val)
{
    // Declare a union to mask the double as an int.
    union AtomicAddUnion
    {
	unsigned long long int d_i;
	double d_d;
	__device__ AtomicAddUnion() {};
    };

    // Create unions for the swap.
    AtomicAddUnion old_value, new_value, next_value;

    // Perform the sum.
    old_value.d_d = *address;
    do 
    {
	next_value.d_i = old_value.d_i;
	new_value.d_d = next_value.d_d + val;
	old_value.d_i = atomicCAS( (unsigned long long int*) address, next_value.d_i , new_value.d_i );
    } while ( next_value.d_i != old_value.d_i );

    return old_value.d_d ;
}
#endif // __NVCC__

//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_Atomic_Add_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Atomic_Add.cuh
//---------------------------------------------------------------------------//
