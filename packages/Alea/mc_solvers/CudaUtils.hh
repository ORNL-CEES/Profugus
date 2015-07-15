//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils.hh
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_CudaUtils_hh
#define mc_solver_CudaUtils_hh

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#endif

namespace alea
{

#ifndef USE_LDG
#define USE_LDG 1
#endif

#ifndef STRUCT_MATRIX
#define STRUCT_MATRIX 1
#endif 

//#ifdef __CUDACC__
//#if STRUCT_MATRIX
struct device_row_data{
	double H;	
	double P;
        double W;
	int inds;
};
//#endif
//#endif

// lower_bound implementation that can be called from device
__device__ inline const double * lower_bound(const double * first,
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

#if STRUCT_MATRIX
// lower_bound implementation that can be called from device
__device__ inline const device_row_data * lower_bound(const device_row_data* first,
                  const device_row_data* last, 
                  double val)
{
       const device_row_data* it;     
       int count, step;
       count = last - first;      
       while( count > 0 )
       {
            step = count/2;
            it = first + step;
#if USE_LDG
            if( __ldg( &(*it).P ) < val ) //Modified by Max
#else
            if( *it.P < val )
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
#endif

// atomicAdd, not provided by Cuda for doubles
__device__ inline double atomicAdd(double* address, double val)
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

__device__ inline void getNewState(int &state, double &wt,
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

#if STRUCT_MATRIX
__device__ inline void getNewState(int &state, double &wt,
              device_row_data* data,
              const int    * const offsets,
              curandState   *rng_state )
{
    // Generate random number
    double rand = curand_uniform_double(rng_state);

    // Sample cdf to get new state 
    auto beg_row = &data[offsets[state]];
    auto end_row = &data[offsets[state+1]];

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
    auto index = elem - data;
#if USE_LDG
    state  =  __ldg( &(data[index].inds) ); //modified by Max
    wt    *=  __ldg( &(data[index].W) ); //modified by Max
#else
    state = data[index].inds;
    wt *= data[index].W;
#endif
}
#endif


__device__ inline void getNewState2(int &state, double &wt,
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
__global__ inline void initialize_rng(curandState *state, int seed, int offset)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed,tid,offset,&state[tid]);

}


__global__ inline void initialize_rng2(curandState *state, int*seed, int offset)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed[tid], 0, offset, &state[tid]);
}



        
}


#endif // mc_solver_CudaUtils_hh
