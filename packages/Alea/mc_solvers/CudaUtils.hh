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

struct device_row_data{
	double H;	
	double P;
        double W;
	int inds;
};


enum class SEED_TYPE{
	SAME,
        DIFF,
        RAND
}; 


class StandardAccess{
public:
      __device__ static inline double get(const double * d){
      	return *d;
      };
 
      __device__ static inline int get(const int * d){
        return *d;
      }
};


class LDGAccess{
public:
      __device__ static inline double get(const double *d){
      	return __ldg(d);
      };
      __device__ static inline int get(const int *d){
        return __ldg(d);
      }
};


// lower_bound implementation that can be called from device
template<class MemoryAccess>
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
        if( MemoryAccess::get(it) < val )
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

// lower_bound implementation that can be called from device
template<class MemoryAccess>
__device__ inline const device_row_data * lower_bound(const device_row_data* first,
                  const device_row_data* last, 
                  double val)
{
       const device_row_data* it;     
       int count, step;
       count = (last - first);      
       while( count > 0 )
       {
            step = count/2;
            it = first + step;
            if( MemoryAccess::get(&it->P) < val ) //Modified by Max
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

template <class MemoryAccess>
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
    auto elem = lower_bound<MemoryAccess>(beg_row,end_row,rand);
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
    state  =  MemoryAccess::get(&inds[index]); //modified by Max
    wt    *=  MemoryAccess::get(&W[index]); //modified by Max
}

template<class MemoryAccess>
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

    auto elem = lower_bound<MemoryAccess>(beg_row,end_row,rand);
    //auto elem = thrust::lower_bound( thrust::seq, beg_row, end_row, rand);

    if( elem == end_row )
    {
        // Invalidate all row data
        state = -1;
        wt = 0.0;
        return;
    }

    // Modify weight and update state
    state  =  MemoryAccess::get( &(elem->inds) ); //modified by Max
    wt    *=  MemoryAccess::get( &(elem->W) ); //modified by Max
}

template<class MemoryAccess>
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
    auto elem = lower_bound<MemoryAccess>(beg_row,end_row,rand);
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
    state  =  MemoryAccess::get(&inds[index]); //modified by Max
    wt    *=  MemoryAccess::get(&W[index]); //modified by Max
}


//---------------------------------------------------------------------------//
/*!
 * \brief Initialize Cuda RNG
 */
//---------------------------------------------------------------------------//
__global__ inline void initialize_rng(curandState *state, int seed, int offset, SEED_TYPE seed_type)
{

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed,tid,offset,&state[tid]);

}


__global__ inline void initialize_rng(curandState *state, int* seed, int offset, SEED_TYPE seed_type)
{

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed[tid], 0, offset, &state[tid]);

}
        
}


#endif // mc_solver_CudaUtils_hh
