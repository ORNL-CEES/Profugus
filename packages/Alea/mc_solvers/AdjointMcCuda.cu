//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcCuda.cu
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

#include "AdjointMcCuda.hh"
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"

#include "CudaUtils.hh"

namespace alea
{

#ifndef USE_LDG
#define USE_LDG 0
#endif

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize history into new state
 */
//---------------------------------------------------------------------------//
__device__ void initializeHistory(int &state, double &wt, int N,
        const double * const start_cdf,
        const double * const start_wt,
              curandState   *rng_state)
{
    // Generate random number
    double rand = curand_uniform_double(rng_state);

    // Sample cdf to get new state
    auto elem = lower_bound(start_cdf,start_cdf+N,rand);

    if( elem == &start_cdf[N-1]+1 )
    {
        state = -1;
        wt    = 0.0;
        return;
    }

    // Get weight and update state
    state = elem-start_cdf;
#if USE_LDG
    wt    = __ldg(&start_wt[state]); //modified by Max
#else
    wt = start_wt[state];
#endif
}


__device__ void initializeHistory2(int &state, double &wt, int N,
        const double * const start_cdf,
        const double * const start_wt,
              double &rand)
{

    // Sample cdf to get new state
    auto elem = lower_bound(start_cdf,start_cdf+N,rand);

    if( elem == &start_cdf[N-1]+1 )
    {
        state = -1;
        wt    = 0.0;
        return;
    }

    // Get weight and update state
    state = elem-start_cdf;
#if USE_LDG
    wt    = __ldg(&start_wt[state]); //modified by Max
#else
    wt = start_wt[state];
#endif
}



//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
__device__ void tallyContribution(int state, double wt,
              double * const x, 
        const double * const H,
        const int    * const inds,
        const int    * const offsets,
              bool           expected_value)
{
    if( expected_value )
    {
        int row_begin = offsets[state];
        int row_end   = offsets[state+1];

        // For expected value estimator, loop over current row and add
        // contributions corresponding to each element
        for( int i=row_begin; i<row_end; ++i )
        {
#if USE_LDG
            atomicAdd(x+inds[i],wt* ( __ldg(&H[i]) ) );//modified by Max
#else
            atomicAdd(x+inds[i],wt*H[i]);//modified by Max
#endif

        }
    }
    else
    {
        // Collision estimator just adds weight
        atomicAdd(x+state,wt);
    }
}


//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
__global__ void adjoint_run_monte_carlo(int N, int history_length, double wt_cutoff,
        int batch_size,
        bool expected_value,
        const double * const start_cdf,
        const double * const start_wt,
        const double * const H,
        const double * const P,
        const double * const W,
        const int    * const inds,
        const int    * const offsets,
        const double * const coeffs,
              double * const x,
              curandState   *rng_state)
{
    int state = -1;
    double wt = 0.0;

    // Store rng state locally
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    curandState local_state = rng_state[tid];
 
    extern __shared__ double steps[];
 
    for (int i = 0; i<batch_size; ++i)
        steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);



    // Get initial state for this history by sampling from start_cdf
    initializeHistory(state,wt,N,start_cdf,start_wt,&local_state);

    //initializeHistory2(state,wt,N,start_cdf,start_wt,steps[threadIdx.x]);

    //printf("Starting history in state %i with weight %7.3f\n",state,wt);
    if( state == -1 )
    {
        rng_state[tid] = local_state;
        return;
    }
    double init_wt = wt;

    // With expected value estimator we start on stage 1 because
    // zeroth order term is added explicitly at the end
    int stage = expected_value ? 1 : 0;

    // Perform initial tally
    tallyContribution(state,coeffs[stage]*wt,x,H,inds,offsets,
        expected_value);

    int count_batch = 0;

    for( ; stage<=history_length; ++stage )
    {
        if (count_batch == batch_size)
        {

          //__syncthreads();
         count_batch = 0;
         for (int i = 0; i<batch_size; ++i)
            steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);
        }

        // Move to new state
        //getNewState(state,wt,P,W,inds,offsets,&local_state);
        //printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

        getNewState2(state,wt,P,W,inds,offsets,steps[threadIdx.x + count_batch * blockDim.x]);

        if( state == -1 )
            break;

        // Tally
        tallyContribution(state,coeffs[stage]*wt,x,H,inds,offsets,
            expected_value);

        count_batch++;

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
AdjointMcCuda::AdjointMcCuda(
        Teuchos::RCP<const MC_Data> mc_data,
        const const_scalar_view     coeffs,
        Teuchos::RCP<Teuchos::ParameterList> pl)

  : d_N(mc_data->getIterationMatrix()->getGlobalNumRows())
{
    // Get parameters
    d_num_histories      = pl->get("num_histories",1000);
    d_max_history_length = coeffs.size()-1;
    d_weight_cutoff      = pl->get("weight_cutoff",0.0);

    // Determine type of tally
    std::string estimator = pl->get<std::string>("estimator",
                                                 "expected_value");
    VALIDATE(estimator == "collision" ||
             estimator == "expected_value",
             "Only collision and expected_value estimators are available.");
    d_use_expected_value = (estimator == "expected_value");

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
void AdjointMcCuda::solve(const MV &b, MV &x)
{
    // Containers to hold starting cdf and wt arrays
    thrust::device_vector<double> start_cdf, start_wt;

    // Build initial probability and weight distributions
    Teuchos::ArrayRCP<const double> b_data = b.getData(0);
    build_initial_distribution(b_data,start_cdf,start_wt);

    // Get pointers for kernel
    const double * const start_cdf_ptr =
        thrust::raw_pointer_cast(start_cdf.data());
    const double * const start_wt_ptr  =
        thrust::raw_pointer_cast(start_wt.data());

    const double * const H       = thrust::raw_pointer_cast(d_H.data());
    const double * const P       = thrust::raw_pointer_cast(d_P.data());
    const double * const W       = thrust::raw_pointer_cast(d_W.data());
    const int    * const inds    = thrust::raw_pointer_cast(d_inds.data());
    const int    * const offsets = thrust::raw_pointer_cast(d_offsets.data());
    const double * const coeffs  = thrust::raw_pointer_cast(d_coeffs.data());

    // Create vector for state
    thrust::device_vector<double> x_vec(d_N);
    double * const x_ptr = thrust::raw_pointer_cast(x_vec.data());

    int block_size = std::min(256,d_num_histories);
    int num_blocks = d_num_histories / block_size;

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

    adjoint_run_monte_carlo<<< num_blocks,block_size, sizeof(double) * block_size * batch_size >>>(d_N,d_max_history_length,
        d_weight_cutoff, batch_size, d_use_expected_value,
        start_cdf_ptr,start_wt_ptr,H,P,W,inds,offsets,coeffs,x_ptr,rng_states);

    // Check for errors in kernel launch
    e = cudaGetLastError();
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to execute MC kernel");

    // Scale by history count
    for( auto itr= x_vec.begin(); itr != x_vec.end(); ++itr )
        *itr /= static_cast<double>(num_blocks*block_size);

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
