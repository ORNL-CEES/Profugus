//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils.cu
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//


#include "CudaUtils.hh"

#include <iterator>
#include <cmath>
#include <curand_kernel.h>
#include <curand.h>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/binary_search.h>
#include <thrust/generate.h>
#include <thrust/random.h>
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"


namespace alea
{


#ifndef USE_LDG
#define USE_LDG 0
#endif


// lower_bound implementation that can be called from device
__device__ const double * lower_bound(const double * first,
        const double * last,
        double   val)
{
    const double * it;
    int count, step;
    count = last - first;
    while( count > 0 )
    {
        step = count / 2;
        it = first+step;
#if USE_LDG
        if( __ldg( &(*it) ) < val ) //Modified by Max
#else
        if ( *it<val )
#endif
        {
            first = ++it;
            count -= step+1;
        }
        else
        {
            count = step;
        }
    }
    return first;
}

// atomicAdd, not provided by Cuda for doubles
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
        (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                __double_as_longlong(val +
                    __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
__device__ void getNewState(int &state, double &wt,
        const double * const P,
        const double * const W,
        const int    * const inds,
        const int    * const offsets,
              curandState   *rng_state )
{
    // Generate random number
    double rand = curand_uniform_double(rng_state);

    // Sample cdf to get new state
    auto beg_row = P + offsets[state];
    auto end_row = P + offsets[state+1];
    auto elem = lower_bound(beg_row,end_row,rand);
    //auto elem = thrust::lower_bound( thrust::seq, beg_row, end_row, rand);

    if( elem == end_row )
    {
        // Invalidate all row data
        state = -1;
        wt = 0.0;
        return;
    }

    // Modify weight and update state
    auto index = elem - P;
#if USE_LDG
    state  =  __ldg(&inds[index]); //modified by Max
    wt    *=  __ldg(&W[index]); //modified by Max
#else
    state = inds[index];
    wt *= W[index];
#endif
}


__device__ void getNewState2(int &state, double &wt,
        const double * const P,
        const double * const W,
        const int    * const inds,
        const int    * const offsets,
              double   &rand )
{

    // Sample cdf to get new state
    auto beg_row = P + offsets[state];
    auto end_row = P + offsets[state+1];
    auto elem = lower_bound(beg_row,end_row,rand);
    //auto elem = thrust::lower_bound( thrust::seq, beg_row, end_row, rand);

    if( elem == end_row )
    {
        // Invalidate all row data
        state = -1;
        wt = 0.0;
        return;
    }

    // Modify weight and update state
    auto index = elem - P;
#if USE_LDG
    state  =  __ldg(&inds[index]); //modified by Max
    wt    *=  __ldg(&W[index]); //modified by Max
#else
    state = inds[index];
    wt *= W[index];
#endif
}


//---------------------------------------------------------------------------//
/*!
 * \brief Initialize Cuda RNG
 */
//---------------------------------------------------------------------------//
__global__ void initialize_rng(curandState *state, int seed, int offset)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed,tid,offset,&state[tid]);

}


__global__ void initialize_rng2(curandState *state, int*seed, int offset)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed[tid], 0, offset, &state[tid]);
}



}
