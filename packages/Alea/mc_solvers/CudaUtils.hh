//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils.hh
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_CudaUtils_hh
#define mc_solver_CudaUtils_hh

#include <cmath>
#include <ctime>

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
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



class OnTheFly
{

public: 
       OnTheFly(curandState*, unsigned int, unsigned int){};
        __device__ inline double get(curandState* rng_state)
       {
              double rand = curand_uniform_double(rng_state);	
	      return rand;
       };
};


__global__ inline void initial_state(curandState* state, double* ss)
{
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
        ss[tid] = curand_uniform_double(&state[tid]);
}



class Precomputed
{

private:
        bool computed = false;
        double* starting_states;
public:
        inline Precomputed(curandState*, unsigned int, unsigned int);
	__device__ inline double get(curandState*)
        { 
          if( computed == false )
            return -1;
          unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x; 
	  return starting_states[tid]; 
        };
};

inline Precomputed::Precomputed(curandState* state, 
	unsigned int num_blocks, unsigned int block_size):computed(true)
{ 
	thrust::device_vector< double > device_ss( num_blocks * block_size );
	thrust::device_ptr< double > dev_ptr = device_ss.data();
	starting_states = thrust::raw_pointer_cast(dev_ptr);
	initial_state<<<num_blocks, block_size>>>(state, starting_states);
	thrust::sort( device_ss.begin(), device_ss.end() );
}

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

class SameSeed{

    public:
          inline SameSeed(unsigned int);
          __device__ inline unsigned int getSeed(){return seed;};
          __device__ inline unsigned int getSequence(){unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x; return tid;};

    private:
          unsigned int seed;
};

inline SameSeed::SameSeed(unsigned int index):seed(index){};


class DifferentSeed{

    public:
          inline DifferentSeed(unsigned int, unsigned int);
          __device__ inline unsigned int getSeed();
          __device__ inline unsigned int getSequence(){return 0;};

    private:
          unsigned int* seed_ptr;
          thrust::device_vector<unsigned int> seeds;
};

inline DifferentSeed::DifferentSeed(unsigned int N, unsigned int index)
{
	seeds.resize(N);
        thrust::sequence( seeds.begin(), seeds.end(), index );
        seed_ptr = thrust::raw_pointer_cast(seeds.data());
}

__device__ inline unsigned int DifferentSeed::getSeed()
{
        unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	return seed_ptr[tid];
}



class RandomSeed{

    public:
          inline RandomSeed(unsigned int);
          __device__ inline unsigned int getSeed();
          __device__ inline unsigned int getSequence(){return 0;};

    private:
          unsigned int* seed_ptr;
          thrust::device_vector<unsigned int> dev_seeds;
};

inline RandomSeed::RandomSeed(unsigned int N)
{
	dev_seeds.resize(N);
        
        thrust::host_vector<int> host_seeds( N );
        std::srand(std::time(0));
    	thrust::generate(host_seeds.begin(), host_seeds.end(), std::rand);
        dev_seeds=host_seeds;
        
        seed_ptr = thrust::raw_pointer_cast(dev_seeds.data());        
}

__device__ inline unsigned int RandomSeed::getSeed()
{
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
	return seed_ptr[tid];
}

template<class seed_type>
__global__ inline void initialize_rng(curandState *state, int offset, seed_type seed)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    curand_init(seed.getSeed(), seed.getSequence(), offset, &state[tid]);	
}

 
}


#endif // mc_solver_CudaUtils_hh
