//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcCuda.cu
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <curand_kernel.h>
#include <curand.h>
#include <cuda_runtime_api.h>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/binary_search.h>
#include <thrust/generate.h>
#include <thrust/random.h>
#include <thrust/sort.h>

#include "AdjointMcCuda.hh"
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"

namespace alea
{

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 256
#endif 

#ifndef BATCH_SIZE
#define BATCH_SIZE 5
#endif

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize history into new state
 */
//---------------------------------------------------------------------------//

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


__global__ void initial_state(curandState* state, double* ss)
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
        Precomputed(curandState*, unsigned int, unsigned int);
	__device__ inline double get(curandState*)
        { 
          if( computed == false )
            return;
          unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x; 
	  return starting_states[tid]; 
        };
};

Precomputed::Precomputed(curandState* state, 
	unsigned int num_blocks, unsigned int block_size):computed(true)
{ 
	thrust::device_vector< double > device_ss(num_blocks * block_size);
	thrust::device_ptr< double > dev_ptr = device_ss.data();
	starting_states = thrust::raw_pointer_cast(dev_ptr);
	initial_state<<<num_blocks, block_size>>>(state, starting_states);
	thrust::sort( device_ss.begin(), device_ss.end() );
}


template<class MemoryAccess, class InitializePolicy>
__device__ void initializeHistory(int &state, double &wt, int N,
        const double * const start_cdf,
        const double * const start_wt,
              curandState   *rng_state,
              InitializePolicy   initialize)
{
    // Generate random number
    double rand = initialize.get( rng_state );

    // Sample cdf to get new state
    auto elem = lower_bound<MemoryAccess>(start_cdf,start_cdf+N,rand);

    if( elem == &start_cdf[N-1]+1 )
    {
        state = -1;
        wt    = 0.0;
        return;
    }

    // Get weight and update state
    state = elem-start_cdf;
    wt    = MemoryAccess::get(&start_wt[state]); //modified by Max
}

template<class MemoryAccess>
__device__ void initializeHistory2(int &state, double &wt, int N,
        const double * const start_cdf,
        const double * const start_wt,
              double &rand)
{

    // Sample cdf to get new state
    auto elem = lower_bound<MemoryAccess>(start_cdf,start_cdf+N,rand);

    if( elem == &start_cdf[N-1]+1 )
    {
        state = -1;
        wt    = 0.0;
        return;
    }

    // Get weight and update state
    state = elem-start_cdf;
    wt    = MemoryAccess::get(&start_wt[state]); //modified by Max
}



//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
template<class MemoryAccess>
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
            atomicAdd(x+inds[i],wt* ( MemoryAccess::get(&H[i]) ) );//modified by Max

    }
    else
    {
        // Collision estimator just adds weight
        atomicAdd(x+state,wt);
    }
}

template<class MemoryAccess>
__device__ void tallyContribution(int state, double wt,
              double * const x, 
              const device_row_data * data,
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
            atomicAdd(x+data[i].inds,
             wt* ( MemoryAccess::get(&(data[i].H)) ));//modified by Max
    }
    else
    {
        // Collision estimator just adds weight
        atomicAdd(x+state,wt);
    }
}

template<class MemoryAccess, class InitializePolicy>
__global__ void run_adjoint_monte_carlo(int N, int history_length, double wt_cutoff,
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
        curandState  * rng_state, 
        InitializePolicy  initialize    )
{
    int state = -1;
    double wt = 0.0;

    // Store rng state locally
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    curandState local_state = rng_state[tid];
 
    //__shared__ double steps[BLOCK_SIZE * BATCH_SIZE];
 
    /*for (int i = 0; i<BATCH_SIZE; ++i)
        steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);
    */
    // Get initial state for this history by sampling from start_cdf
    initializeHistory< MemoryAccess >(state,wt,N,start_cdf,start_wt,&local_state, initialize);

    //initializeHistory2< MemoryAccess >(state,wt,N,start_cdf,start_wt,steps[threadIdx.x]);

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
    tallyContribution< MemoryAccess >(state,coeffs[stage]*wt,x,H,inds,offsets,
        expected_value);

    //int count_batch = 1;

    for( ; stage<=history_length; ++stage )
    {
    /*    if (count_batch == BATCH_SIZE)
        {

          //__syncthreads();
         count_batch = 0;
         for (int i = 0; i<BATCH_SIZE; ++i)
            steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);
        }*/

        // Move to new state
        getNewState< MemoryAccess >(state,wt,P,W,inds,offsets,&local_state);
        //printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

	//getNewState2< MemoryAccess >(state,wt,P,W,inds,offsets,steps[threadIdx.x + count_batch * blockDim.x]);

        if( state == -1 )
            break;

        // Tally
        tallyContribution< MemoryAccess >(state,coeffs[stage]*wt,x,H,inds,offsets,
            expected_value);

       // count_batch++;

        // Check weight cutoff
        if( std::abs(wt/init_wt) < wt_cutoff )
            break;
   
    }

    // Store rng state back to global
    rng_state[tid] = local_state;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
template<class MemoryAccess, class InitializePolicy>
__global__ void run_adjoint_monte_carlo(int N, int history_length, double wt_cutoff,
        bool expected_value,
        const double    * const start_cdf,
        const double    * const start_wt,
        device_row_data * data,
        const int       * const offsets,
        const double    * const coeffs,
              double    * const x,
        curandState     * rng_state,
        InitializePolicy     initialize  )
{
    int state = -1;
    double wt = 0.0;

    // Store rng state locally
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    curandState local_state = rng_state[tid];
 
    // Get initial state for this history by sampling from start_cdf
    initializeHistory< MemoryAccess >(state,wt,N,start_cdf,start_wt,&local_state, initialize);

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
    tallyContribution< MemoryAccess >(state,coeffs[stage]*wt,x,data,offsets,
        expected_value);

    for( ; stage<=history_length; ++stage )
    {

        // Move to new state
        getNewState< MemoryAccess >(state,wt,data,offsets,&local_state);
        //printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

        if( state == -1 )
            break;

        // Tally
        tallyContribution< MemoryAccess >(state,coeffs[stage]*wt,x,data,offsets,
            expected_value);

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
    d_struct             = pl->get("struct_matrix", 0);
    d_use_ldg            = pl->get("use_ldg", 0);
    d_precompute_states  = pl->get("precompute_states", 0);
    std::string seed_type= pl->get("seed_type", std::string("same") );

    VALIDATE( d_struct==0 || d_struct==1, 
            "Value for the flag to manage matrix data not valid" );
            
    VALIDATE( d_use_ldg==0 || d_use_ldg==1, 
            "Value for the texture memory handling not valid" );   
            
    VALIDATE( d_precompute_states==0 || d_precompute_states==1, 
            "Value for the policy of computing states not valid" );                  

    VALIDATE( seed_type == std::string("same") 
              || seed_type == std::string("different")
              || seed_type == std::string("random"), 
              "Type of seed selected is not valid" );	

    if( seed_type.c_str()==std::string("same") )   
    	d_seed_type = SEED_TYPE::SAME;
    else if( seed_type.c_str()==std::string("different") )
        d_seed_type = SEED_TYPE::DIFF;
    else if( seed_type.c_str()==std::string("random") ) 
        d_seed_type = SEED_TYPE::RAND;

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



template <class InitializePolicy>
void AdjointMcCuda::launch_monte_carlo(int d_N,int num_blocks,
     int d_max_history_length, double d_weight_cutoff, bool d_use_expected_value,
     const double * const start_cdf_ptr,const double * const start_wt_ptr,
     const double * H,
     const double * P,
     const double * W,
     const int * inds,
     device_row_data * data_ptr,
     const int * const offsets,
     const double * const coeffs,
     double * x_ptr,
     curandState * rng_states)
{

    InitializePolicy initialize( rng_states, num_blocks,  BLOCK_SIZE );

    if( d_struct==0 )
    {
        if( d_use_ldg==0 )
        {
            run_adjoint_monte_carlo<StandardAccess><<< num_blocks,BLOCK_SIZE >>>
		   (d_N, d_max_history_length, d_weight_cutoff, 
		   d_use_expected_value,
                   start_cdf_ptr,start_wt_ptr,H,P,W,
                   inds,offsets,coeffs,x_ptr,rng_states, initialize);
        }
        else
        {
            run_adjoint_monte_carlo<LDGAccess><<< num_blocks,BLOCK_SIZE >>>
                   (d_N,d_max_history_length, d_weight_cutoff,
                   d_use_expected_value,
                   start_cdf_ptr,start_wt_ptr,H,P,W,
                   inds,offsets,coeffs,x_ptr,rng_states, initialize);

        }        
    }
    else
    {
        if( d_use_ldg==0 )
        {
            run_adjoint_monte_carlo<StandardAccess><<< num_blocks,BLOCK_SIZE >>>(d_N,
                    d_max_history_length, d_weight_cutoff, d_use_expected_value,
                    start_cdf_ptr,start_wt_ptr,data_ptr,
                    offsets,coeffs,x_ptr,rng_states, initialize);
        }
        else{
            run_adjoint_monte_carlo<LDGAccess><<< num_blocks,BLOCK_SIZE >>>(d_N,
                    d_max_history_length, d_weight_cutoff, d_use_expected_value,
                    start_cdf_ptr,start_wt_ptr,data_ptr,
                    offsets,coeffs,x_ptr,rng_states, initialize);

        }	                   
    }          
 

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

    const double * H;
    const double * P;
    const double * W;
    const int    * inds; 
    device_row_data * data_ptr;

    if( d_struct==0 )
    {
    	H       = thrust::raw_pointer_cast(d_H.data());
    	P       = thrust::raw_pointer_cast(d_P.data());
    	W       = thrust::raw_pointer_cast(d_W.data());
    	inds    = thrust::raw_pointer_cast(d_inds.data());
    }
    else
    	 data_ptr = thrust::raw_pointer_cast(mat_data.data());   

    const int    * const offsets = thrust::raw_pointer_cast(d_offsets.data());
    const double * const coeffs  = thrust::raw_pointer_cast(d_coeffs.data());

    // Create vector for state
    thrust::device_vector<double> x_vec(d_N);
    double * const x_ptr = thrust::raw_pointer_cast(x_vec.data());

    //int block_size = std::min(256,d_num_histories);
    VALIDATE( BLOCK_SIZE <= d_num_histories, 
          "Number of histories is smaller than the block size" );

    int num_blocks = d_num_histories / BLOCK_SIZE;

    curandState *rng_states;
    cudaError e = cudaMalloc( (void **)&rng_states,
        BLOCK_SIZE * num_blocks * sizeof(curandState) );

    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to allocate memory");

    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    if( d_seed_type==SEED_TYPE::SAME )
    {
        std::cout<<"Same seed instantiated for all the threads"<<std::endl;

    	// Initialize RNG
    	SameSeed seed(d_rng_seed);
	initialize_rng<SameSeed><<<num_blocks,BLOCK_SIZE>>>(rng_states,
    	    d_num_curand_calls, seed);
    }
    else if ( d_seed_type==SEED_TYPE::DIFF )
    {
        std::cout<<"Different adjacent seeds instantiated"<<std::endl;

 	DifferentSeed seed( BLOCK_SIZE*num_blocks, d_rng_seed );
    	
    	initialize_rng<DifferentSeed><<<num_blocks, BLOCK_SIZE>>>(rng_states,  
            d_num_curand_calls, seed);
    }
    else if ( d_seed_type==SEED_TYPE::RAND )
    {
        std::cout<<"Different random seeds instantiated from 0 to "<<
         RAND_MAX<<std::endl;

    	RandomSeed seed( BLOCK_SIZE*num_blocks);
 
    	initialize_rng<RandomSeed><<<num_blocks, BLOCK_SIZE>>>(rng_states,
            d_num_curand_calls, seed);
    }

    // Check for errors in kernel launch
    e = cudaGetLastError();
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to initialize RNG");
    d_num_curand_calls++;

    if( d_precompute_states == 0 )
        launch_monte_carlo< OnTheFly >(d_N,num_blocks,
              d_max_history_length, d_weight_cutoff, d_use_expected_value,
              start_cdf_ptr,start_wt_ptr,H,P,W,
              inds,data_ptr,offsets,coeffs,x_ptr,rng_states);
    else
        launch_monte_carlo< Precomputed >(d_N,num_blocks,
              d_max_history_length, d_weight_cutoff, d_use_expected_value,
              start_cdf_ptr,start_wt_ptr,H,P,W,
              inds,data_ptr,offsets,coeffs,x_ptr,rng_states);

    // Scale by history count
    for( auto itr= x_vec.begin(); itr != x_vec.end(); ++itr )
        *itr /= static_cast<double>( num_blocks * BLOCK_SIZE );

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
        d_offsets.resize(d_N+1);

	if(d_struct == 0)
	{
   		d_nnz = H->getNodeNumEntries();
    		d_H.resize(d_nnz);
    		d_P.resize(d_nnz);
    		d_W.resize(d_nnz);
    		d_inds.resize(d_nnz);
	
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

	else
	{
                d_nnz = H->getNodeNumEntries();
    		Teuchos::ArrayView<const double> pval_row;
    		Teuchos::ArrayView<const double> hval_row;
    		Teuchos::ArrayView<const double> wval_row;
    		Teuchos::ArrayView<const int>    ind_row;
    // This loop should perhaps be rewritten, right now a separate call
    // to cudaMemcpy is performed for each row of each matrix
    // It might be more efficient to create a single vector on the CPU
    // and do a single copy to device?
    		d_offsets[0] = 0;

    		thrust::host_vector< device_row_data > data_host( d_nnz );
    		mat_data.resize( d_nnz );

    		int count = 0;
    		for( int i=0; i<d_N; ++i )
    		{
        		// Extract row i of matrix
        		H->getLocalRowView(i,ind_row,hval_row);
        		P->getLocalRowView(i,ind_row,pval_row);
        		W->getLocalRowView(i,ind_row,wval_row); 
     
        		for( int j = 0; j < ind_row.size(); ++j )
        		{
                    CHECK( count < data_host.size() );
                    CHECK( j < hval_row.size() );
                    CHECK( j < pval_row.size() );
                    CHECK( j < wval_row.size() );
                    CHECK( j < ind_row.size() );
                    data_host[count].H   = hval_row[j];            
                    data_host[count].P   = pval_row[j];
                    data_host[count].W   = wval_row[j];
                    data_host[count].inds = ind_row[j];
                    count++;
        		}
        
        		d_offsets[i+1] = d_offsets[i] + ind_row.size();
    		}

    		CHECK( count == d_nnz );

    		thrust::copy( data_host.begin(), data_host.end(), mat_data.begin() );
            //mat_data = data_host;

    		// Copy coefficients into device vector
    		const_scalar_view::HostMirror coeffs_host = Kokkos::create_mirror_view(coeffs);
    		Kokkos::deep_copy(coeffs_host,coeffs);
    		d_coeffs.resize(coeffs.size());
    		thrust::copy(coeffs_host.ptr_on_device(),
                 	coeffs_host.ptr_on_device()+coeffs_host.size(),
                 	d_coeffs.begin());

       }
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
