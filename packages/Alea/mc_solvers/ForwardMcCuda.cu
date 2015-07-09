//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ForwardMcCuda.cu
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

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

#include "ForwardMcCuda.hh"
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"

namespace alea
{

#define USE_LDG 0

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

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize history into new state
 */
//---------------------------------------------------------------------------//
__device__ void initializeHistory(int &state, const double * const b, double &wt)
{
  wt = b[state];
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
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
__device__ void tallyContribution(int state, double wt, double * const x)
{
        // Collision estimator just adds weight
        atomicAdd(x+state,wt);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
__global__ void run_monte_carlo(int N, int history_length, double wt_cutoff,
        int entry_histories, 
        int batch_size,
        const double * const H,
        const double * const P,
        const double * const W,
        const int    * const inds,
        const int    * const offsets,
        const double * const coeffs,
              double * const x,
        const double * const rhs, 
              curandState   *rng_state)
{
    int state = -1;
    double wt = 0.0;

    // Store rng state locally
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int entry = tid % entry_histories;
    state = entry;
    curandState local_state = rng_state[tid];
 
/*    extern __shared__ double steps[];
 
    for (int i = 0; i<batch_size; ++i)
        steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);

*/

    // Get initial state for this history by sampling from start_cdf
    initializeHistory(state, rhs, wt);

    //initializeHistory2(state,wt,N,start_cdf,start_wt,steps[threadIdx.x]);

    //printf("Starting history in state %i with weight %7.3f\n",state,wt);
    if( state > N )
    {
        rng_state[tid] = local_state;
        return;
    }
    double init_wt = wt;

    // Perform initial tally
    tallyContribution(state,coeffs[stage]*wt,x,H,inds,offsets);

  //  int count_batch = 0;

    for( ; stage<=history_length; ++stage )
    {
/*        if (count_batch == batch_size)
        {

          //__syncthreads();
         count_batch = 0;
         for (int i = 0; i<batch_size; ++i)
            steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);
        }

*/

        // Move to new state
        getNewState(state,wt,P,W,inds,offsets,&local_state);
        //printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

        //getNewState2(state,wt,P,W,inds,offsets,steps[threadIdx.x + count_batch * blockDim.x]);

        if( state == -1 )
            break;

        // Tally
        tallyContribution(entry,coeffs[stage]*wt,x);

        //count_batch++;

        // Check weight cutoff
        if( std::abs(wt/init_wt) < wt_cutoff )
            break;
   
    }

    // Store rng state back to global
    rng_state[tid] = local_state;
}
        
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param P Views into entries of probability matrix
 * \param W Views into entries of weight matrix
 * \param inds Views into nonzeros indices
 * \param offsets Starting indices for each matrix row
 * \param coeffs Polynomial coefficients
 * \param pl Problem parameters
 */
//---------------------------------------------------------------------------//
ForwardMcCuda::ForwardMcCuda(
        Teuchos::RCP<const MC_Data> mc_data,
        const const_scalar_view     coeffs,
        Teuchos::RCP<Teuchos::ParameterList> pl)

  : d_N(mc_data->getIterationMatrix()->getGlobalNumRows())
{
    // Get parameters
    d_num_histories      = pl->get("num_histories",1000);
    d_max_history_length = coeffs.size()-1;
    d_weight_cutoff      = pl->get("weight_cutoff",0.0);

    VALIDATE(estimator == "collision",
             "Only collision estimator is available.");

    // Should we print anything to screen
    std::string verb = profugus::to_lower(pl->get("verbosity","low"));
    if( verb == "none" )
        d_verbosity = NONE;
    else if( verb == "low" )
        d_verbosity = LOW;
    else if( verb == "high" )
        d_verbosity = HIGH;

    prepareDeviceData(mc_data,coeffs);

    d_num_curand_calls = 0;
    d_rng_seed = pl->get<int>("rng_seed",1234);
}

//---------------------------------------------------------------------------//
// Solve problem using Monte Carlo
//---------------------------------------------------------------------------//
void ForwardMcCuda::solve(const MV &b, MV &x)
{
    Teuchos::ArrayRCP<const double> b_data = b.getData(0);
    const double * const rhs = b.getRawPtr();

    const double * const H       = thrust::raw_pointer_cast(d_H.data());
    const double * const P       = thrust::raw_pointer_cast(d_P.data());
    const double * const W       = thrust::raw_pointer_cast(d_W.data());
    const int    * const inds    = thrust::raw_pointer_cast(d_inds.data());
    const int    * const offsets = thrust::raw_pointer_cast(d_offsets.data());
    const double * const coeffs  = thrust::raw_pointer_cast(d_coeffs.data());

    // Create vector for state
    thrust::device_vector<double> x_vec(d_N);
    double * const x_ptr = thrust::raw_pointer_cast(x_vec.data());

    int tot_histories = d_num_histories * d_N;

    int block_size = std::min(256,tot_histories);
    int num_blocks = tot_histories / block_size;
    
    VALIDATE( num_blocks > 0, "The size of the problem is too small" );

    curandState *rng_states;
    cudaError e = cudaMalloc((void **)&rng_states,
        block_size*num_blocks*sizeof(curandState));

    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to allocate memory");

    // Initialize RNG
    //initialize_rng<<<num_blocks,block_size>>>(rng_states,d_rng_seed,
    //    d_num_curand_calls);

    thrust::device_vector<int> seeds( block_size*num_blocks);
    thrust::sequence(seeds.begin(), seeds.end(), d_rng_seed);
    int* seed_ptr = thrust::raw_pointer_cast(seeds.data());

    initialize_rng2<<<num_blocks, block_size>>>(rng_states, seed_ptr, 
          d_num_curand_calls);

    // Check for errors in kernel launch
    e = cudaGetLastError();
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to initialize RNG");
    d_num_curand_calls++;

    int batch_size = 5;

    run_monte_carlo<<< num_blocks,block_size, sizeof(double) * block_size * batch_size >>>(d_N,d_max_history_length,
        d_weight_cutoff, d_num_histories, batch_size,
        H,P,W,inds,offsets,coeffs,x_ptr, rhs, rng_states);

    // Check for errors in kernel launch
    e = cudaGetLastError();
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to execute MC kernel");

    // Scale by history count
    for( auto itr= x_vec.begin(); itr != x_vec.end(); ++itr )
        *itr /= static_cast<double>(d_num_histories);

    // Copy data back to host
    {
        Teuchos::ArrayRCP<double> x_data = x.getDataNonConst(0);
        thrust::copy(x_vec.begin(),x_vec.end(),x_data.get());
    }


    // Add rhs for expected value
    if( d_use_expected_value )
       x.update(d_coeffs[0],b,1.0);

    // Free RNG state
    e = cudaFree(rng_states);
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to deallocate memory");
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Extract matrices into ArrayView objects for faster data access
//---------------------------------------------------------------------------//
void AdjointMcCuda::prepareDeviceData(Teuchos::RCP<const MC_Data> mc_data,
        const const_scalar_view coeffs)
{
    Teuchos::RCP<const MATRIX> H = mc_data->getIterationMatrix();
    Teuchos::RCP<const MATRIX> P = mc_data->getProbabilityMatrix();
    Teuchos::RCP<const MATRIX> W = mc_data->getWeightMatrix();

    d_nnz = H->getNodeNumEntries();
    d_H.resize(d_nnz);
    d_P.resize(d_nnz);
    d_W.resize(d_nnz);
    d_inds.resize(d_nnz);
    d_offsets.resize(d_N+1);

    Teuchos::ArrayView<const double> val_row;
    Teuchos::ArrayView<const int>    ind_row;
    auto h_iter   = d_H.begin();
    auto p_iter   = d_P.begin();
    auto w_iter   = d_W.begin();
    auto ind_iter = d_inds.begin();
    // This loop should perhaps be rewritten, right now a separate call
    // to cudaMemcpy is performed for each row of each matrix
    // It might be more efficient to create a single vector on the CPU
    // and do a single copy to device?
    d_offsets[0] = 0;
    for( int i=0; i<d_N; ++i )
    {
        // Extract row i of matrix
        H->getLocalRowView(i,ind_row,val_row);
        thrust::copy(val_row.begin(),val_row.end(),h_iter);
        h_iter += val_row.size();
        P->getLocalRowView(i,ind_row,val_row);
        thrust::copy(val_row.begin(),val_row.end(),p_iter);
        p_iter += val_row.size();
        W->getLocalRowView(i,ind_row,val_row);
        thrust::copy(val_row.begin(),val_row.end(),w_iter);
        w_iter += val_row.size();
        thrust::copy(ind_row.begin(),ind_row.end(),ind_iter);
        ind_iter += ind_row.size();
        d_offsets[i+1] = d_offsets[i] + ind_row.size();
    }
    CHECK( h_iter   == d_H.end() );
    CHECK( p_iter   == d_P.end() );
    CHECK( w_iter   == d_W.end() );
    CHECK( ind_iter == d_inds.end() );

    // Copy coefficients into device vector
    const_scalar_view::HostMirror coeffs_host = Kokkos::create_mirror_view(coeffs);
    Kokkos::deep_copy(coeffs_host,coeffs);
    d_coeffs.resize(coeffs.size());
    thrust::copy(coeffs_host.ptr_on_device(),
                 coeffs_host.ptr_on_device()+coeffs_host.size(),
                 d_coeffs.begin());
}

//---------------------------------------------------------------------------//
// Build initial cdf and weights
//---------------------------------------------------------------------------//
void AdjointMcCuda::build_initial_distribution(
        Teuchos::ArrayRCP<const double> b,
        thrust::device_vector<double>  &cdf,
        thrust::device_vector<double>  &wt) const
{
    thrust::host_vector<double> cdf_host(d_N);
    thrust::host_vector<double> wt_host(d_N);

    // First take absolute value of b
    for( int i=0; i<d_N; ++i )
    {
        cdf_host[i] = std::abs(b[i]);
    }

    // Normalize to get a PDF
    double pdf_sum = std::accumulate(cdf_host.begin(),cdf_host.end(),0.0);
    ENSURE( pdf_sum > 0.0 );
    std::transform(cdf_host.begin(),cdf_host.end(),cdf_host.begin(),
            [pdf_sum](double val){return val/pdf_sum;});

    // Compute weight vector s.t. b = p * wt_host
    std::transform(b.begin(),b.end(),cdf_host.begin(),wt_host.begin(),
            [](double u, double v){return v==0.0 ? 0.0 : u/v;});

    // Convert PDF to CDF
    std::partial_sum(cdf_host.begin(),cdf_host.end(),cdf_host.begin());

    // Copy data to device
    cdf = cdf_host;
    wt  = wt_host;
}

} // namespace alea
