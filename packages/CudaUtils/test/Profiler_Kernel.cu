//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Profiler_Kernel.cu
 * \author Seth R Johnson
 * \date   Tue Jul 09 08:29:11 2013
 * \brief  Profiler_Kernel kernel definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../cuda_utils/Kernel_Header.cuh"

#include "Profiler_Kernel.cuh"

#include "harness/DBC.hh"
#include "../cuda_utils/CudaDBC.hh"

#ifndef NDEBUG
// debug mode is slower
#define PROFILER_INNER_LOOPS  4
#else
#define PROFILER_INNER_LOOPS  16
#endif

namespace cuda
{
//---------------------------------------------------------------------------//
// FUNCTORS
//---------------------------------------------------------------------------//
#define MAD_16(x, y, z) \
    z = x * y + z; z = x * y + z; z = x * y + z; z = x * y + z; \
    z = x * y + z; z = x * y + z; z = x * y + z; z = x * y + z; \
    z = x * y + z; z = x * y + z; z = x * y + z; z = x * y + z; \
    z = x * y + z; z = x * y + z; z = x * y + z; z = x * y + z;

#define MAD_256(x, y, z) \
    MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); \
    MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); \
    MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); \
    MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z); MAD_16(x, y, z);

template<typename T>
struct Profiler_Op
{
    enum { NUM_OPERATIONS = 16 * 16 * 2 }; // 16 * 16 * (1 [add] + 1 [mul])

    __device__ void operator() (T& a, T& b)
    {
        MAD_256(a, b, a); // a = a * b + a
    }
};

#if 0
/*
 * You can override the profiling operation for a specific type by defining the
 * appropriate functor.
 */
template<>
struct Profiler_Op<int>
{
    __device__ void operator() (T& a, T& b)
    {
    }
};
#endif

//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

template<typename Arch_T, typename Float_T, class Operation>
__global__ void operation_test_kernel(
        Float_T*     buffer,
        unsigned int size,
        Operation    op)
{
    typedef Float_T value_type;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = gridDim.x * blockDim.x;

    for (; i < size; i += stride)
    {
        // Values on which to operate
        value_type a = i;
        value_type b = buffer[i];

        for (unsigned int iter = 0; iter != PROFILER_INNER_LOOPS; ++iter)
        {
            // 16 calls to op(a, b)
            op(a, b); op(a, b); op(a, b); op(a, b);
            op(a, b); op(a, b); op(a, b); op(a, b);
            op(a, b); op(a, b); op(a, b); op(a, b);
            op(a, b); op(a, b); op(a, b); op(a, b);
        }
        buffer[i] = a;
    }
}

//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

/*!
 * \brief Do a bunch of operations on the GPU
 *
 * \return number of operations being performed.
 */
template<typename Arch_T, typename Float_T>
__host__ unsigned long int operation_test(
        unsigned int num_blocks,
        unsigned int num_threads,
        Device_Vector<Arch_T, Float_T>& data)
{
    typedef Profiler_Op<Float_T> Operation_t;

    REQUIRE(num_threads * num_blocks > 0);

    // Functor to run on the device
    Operation_t op;

    // Launch the kernel, and raise an assertion if the launch failed
    operation_test_kernel<Arch_T, Float_T, Operation_t>
#ifdef __NVCC__
        <<<num_blocks, num_threads>>>
#endif
        (data.data(), data.size(), op);

    CudaCall(cudaGetLastError());

    // Total number of operations called
    return 1ul
        * Operation_t::NUM_OPERATIONS // Operations per call to op()
        * 16                        // Number of operations in the inner loop
        * PROFILER_INNER_LOOPS      // Number of inner loops
        * data.size();              // Number of elements (blocks*threads)
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATION
//---------------------------------------------------------------------------//
#ifdef __NVCC__
typedef ::cuda::arch::Device Arch_t;
#else
typedef ::cuda::arch::Host Arch_t;
#endif

//---------------------------------------------------------------------------//
// These will automatically instantiate the global kernels and other necessary
// device functions.
//---------------------------------------------------------------------------//
template __host__ unsigned long int operation_test(unsigned int, unsigned int,
        Device_Vector<Arch_t, float>&);
template __host__ unsigned long int operation_test(unsigned int, unsigned int,
        Device_Vector<Arch_t, double>&);
template __host__ unsigned long int operation_test(unsigned int, unsigned int,
        Device_Vector<Arch_t, int>&);

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Profiler_Kernel.cu
//---------------------------------------------------------------------------//
