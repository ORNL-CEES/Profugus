// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Pseudo_Cuda.hh
 * \author Seth R Johnson
 * \date   Tue Aug 13 13:13:37 2013
 * \brief  Plug-in replacement for CUDA to allow GPU code to run on CPU
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * To access the static variables in this file from a CUDA function,
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Pseudo_Cuda_hh
#define CudaUtils_cuda_utils_Pseudo_Cuda_hh

#ifdef __NVCC__
#error "This file should never be compiled with the CUDA compiler!"
#endif

#ifdef __CUDACC__
#error "This is not a cuda file!"
#endif

//! For some of our files (like CudaDBC) we want to change behavior.
#define PSEUDO_CUDA

//---------------------------------------------------------------------------//
// FUNCTION TYPE QUALIFIERS
//---------------------------------------------------------------------------//
#define __device__
#define __global__
#define __host__

//---------------------------------------------------------------------------//
// VARIABLE TYPE QUALIFIERS
//---------------------------------------------------------------------------//

#define __constant__
#define __shared__
#define __restrict__

//---------------------------------------------------------------------------//
// BUILT-IN VECTOR TYPES
//
// Note: we're only including the most common for now; if later implementations
// need other vector types, you can add them.
//---------------------------------------------------------------------------//
struct uint3
{
    unsigned int x, y, z;
};

//! Like uint3, but the constructor initializes non-specified components to 1
struct dim3
{
    unsigned int x, y, z;

    dim3(unsigned int x_=1, unsigned int y_=1, unsigned int z_=1)
        : x(x_), y(y_), z(z_)
    { /* * */ }
};

//! Float with 3 components
struct float3
{
    float x, y, z;
};

//! Double with 3 components
struct double3
{
    double x, y, z;
};

//---------------------------------------------------------------------------//
// BUILT-IN VARIABLES
//
// Note that we use some inline functions and preprocessor definitions to
// emulate the CUDA structures without using global external variables,
// allowing the block sizes to be inlined.
//---------------------------------------------------------------------------//
namespace pseudocuda_detail
{
//---------------------------------------------------------------------------//
//! Return the dimensions (1,1,1)
inline dim3 get_unit_dim3()
{
    return dim3();
}

inline uint3 get_zeroed_uint3()
{
    uint3 result = {0,0,0};
    return result;
}

inline int get_unit_int()
{
    return 1;
}
//---------------------------------------------------------------------------//
}
//---------------------------------------------------------------------------//

//! Dimensions of the grid (1,1,1)
#define gridDim (::pseudocuda_detail::get_unit_dim3())
//! Index of the block in the thread (0,0,0)
#define blockIdx (::pseudocuda_detail::get_zeroed_uint3())
//! Dimensions of the block (1,1,1)
#define blockDim (::pseudocuda_detail::get_unit_dim3())
//! Index of the thread (0,0,0)
#define threadIdx (::pseudocuda_detail::get_zeroed_uint3())
//! Number of threads per warp (1)
#define warpSize (::pseudocuda_detail::get_unit_int())

//---------------------------------------------------------------------------//
// BUILT-IN FUNCTIONS
//---------------------------------------------------------------------------//

inline void __threadfence() { /* * */ }
inline void __threadfence_block() { /* * */ }
inline void __syncthreads() { /* * */ }

template<typename T>
inline T atomicCAS(T* address, T compare, T val)
{
    T old = *address;
    *address = (old == compare ? val : old);
    return old;
}

//---------------------------------------------------------------------------//
#endif // CudaUtils_cuda_utils_Pseudo_Cuda_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Pseudo_Cuda.hh
//---------------------------------------------------------------------------//
