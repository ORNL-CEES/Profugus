//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ForwardMcCuda.cu
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <cmath>
#include <ctime>
#include <curand_kernel.h>
#include <curand.h>
#include <cuda_profiler_api.h>
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

__device__ void tallyContribution2(double wt, double * const x)
{
        atomicAdd(x,wt);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//

template<class MemoryAccess>
__global__ void run_forward_monte_carlo(int N, int history_length, double wt_cutoff,
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
    double wt = 1.0;

    // Store rng state locally
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int entry = tid / entry_histories;

    if(entry >= N)
      return;

    state = entry;
    curandState local_state = rng_state[tid];
 
/*    extern __shared__ double steps[];
 
    for (int i = 0; i<batch_size; ++i)
        steps[threadIdx.x + i*blockDim.x] = curand_uniform_double(&local_state);

*/


    //initializeHistory2(state,wt,N,start_cdf,start_wt,steps[threadIdx.x]);

    //printf("Starting history in state %i with weight %7.3f\n",state,wt);
    if( state >= N )
    {
        rng_state[tid] = local_state;
        return;
    }
    double init_wt = wt;

    int stage = 0;

    // Perform initial tally
    tallyContribution(state,coeffs[stage]*wt*rhs[state],x);

  //  int count_batch = 0;

    for(; stage<=history_length; ++stage )
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
        getNewState<MemoryAccess>(state,wt,P,W,inds,offsets,&local_state);
        //printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

        //getNewState2<MemoryAccess>(state,wt,P,W,inds,offsets,steps[threadIdx.x + count_batch * blockDim.x]);

        if( state == -1 )
            break;

        // Tally
        tallyContribution(entry,coeffs[stage]*wt*rhs[state],x);

        //count_batch++;

        // Check weight cutoff
        if( std::abs(wt/init_wt) < wt_cutoff )
            break;
   
    }

    // Store rng state back to global
    rng_state[tid] = local_state;
}


template<class MemoryAccess>
__global__ void run_forward_monte_carlo2(int N, int history_length, double wt_cutoff,
        int entry_histories, 
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
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

//    extern __shared__ double sol[];

//    sol[threadIdx.x]=0.0;

    if( tid < N )
    {
        int entry=tid;
        
        for(int i=0; i<entry_histories; ++i)
        {
	    int state = entry;
	    double wt = 1.0;

	    curandState local_state = rng_state[tid];
	 
	    //printf("Starting history in state %i with weight %7.3f\n",state,wt);
	    if( state >= N )
	    {
		rng_state[tid] = local_state;
		return;
	    }
	    double init_wt = wt;

	    int stage = 0;

	    // Perform initial tally
	    tallyContribution(state,coeffs[stage]*wt*rhs[state],x);

            //tallyContribution2(coeffs[stage]*wt*rhs[state],&sol[threadIdx.x]);

	    for(; stage<=history_length; ++stage )
	    {
		// Move to new state
		getNewState<MemoryAccess>(state,wt,P,W,inds,offsets,&local_state);
		//printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

		if( state == -1 )
		    break;

		// Tally
		tallyContribution(entry,coeffs[stage]*wt*rhs[state],x);
                //tallyContribution2(coeffs[stage]*wt*rhs[state],&sol[threadIdx.x]);

		// Check weight cutoff
		if( std::abs(wt/init_wt) < wt_cutoff )
		    break;
	   
	    }

	    // Store rng state back to global
	    rng_state[tid] = local_state;
        }
    
        //x[entry]=sol[threadIdx.x];
        //v[entry]=sol[threaddx.x]*sol[threadIdx.x];

    }
}


template<class MemoryAccess>
__global__ void run_forward_monte_carlo3(int N, 
        int history_length, 
        double wt_cutoff,
        int entry_histories, 
        device_row_data* data, 
        const double * const coeffs,
        const int    * const offsets,
              double * const x,
        const double * const rhs, 
              curandState   *rng_state)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if( tid < N )
    {
        int entry=tid;
        
        for(int i=0; i<entry_histories; ++i)
        {
	    int state = entry;
	    double wt = 1.0;

	    curandState local_state = rng_state[tid];
	 
	    //printf("Starting history in state %i with weight %7.3f\n",state,wt);
	    if( state >= N )
	    {
		rng_state[tid] = local_state;
		return;
	    }
	    double init_wt = wt;

	    int stage = 0;

	    // Perform initial tally
	    tallyContribution(state,coeffs[stage]*wt*rhs[state],x);

            //tallyContribution2(coeffs[stage]*wt*rhs[state],&sol[threadIdx.x],&var[threadIdx.x], variance);

	    for(; stage<=history_length; ++stage )
	    {
		// Move to new state
		getNewState<MemoryAccess>(state,wt,data,offsets,&local_state);
		//printf("Stage %i, moving to state %i with new weight of %7.3f\n",stage,state,wt);

		if( state == -1 )
		    break;

		// Tally
		tallyContribution(entry,coeffs[stage]*wt*rhs[state],x);
                //tallyContribution2(coeffs[stage]*wt*rhs[state],&sol[threadIdx.x], &var[threadIdx.x], variance);

		// Check weight cutoff
		if( std::abs(wt/init_wt) < wt_cutoff )
		    break;
	   
	    }

	    // Store rng state back to global
	    rng_state[tid] = local_state;
        }
    
        //x[entry]=sol[threadIdx.x];

    }
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
    d_num_histories        = pl->get("num_histories",1000);
    d_max_history_length   = coeffs.size()-1;
    d_weight_cutoff        = pl->get("weight_cutoff",0.0);
    d_compute_variance     = pl->get("compute_variance", 0);
    d_struct               = pl->get("struct_matrix", 0);
    d_use_ldg              = pl->get("use_ldg", 0);
    d_device_number        = pl->get("device_number", 0);
    d_use_thread_per_entry = pl->get("thread_per_entry",0);
    std::string seed_type  = pl->get("seed_type", std::string("same"));

    VALIDATE( d_struct == 0 || d_struct == 1, 
            "Value for the flag to manage matrix data not valid" );
           
    VALIDATE( d_use_ldg==0 || d_use_ldg==1, 
            "Value for the texture memory handling not valid" );                    
    VALIDATE( d_use_thread_per_entry==0 || d_use_thread_per_entry==1,
            "Value for the task distribution between threads not valid" );

    VALIDATE( seed_type == std::string("same") 
              || seed_type == std::string("different")
              || seed_type == std::string("random"), 
              "Type of seed selected is not valid" );	

    int device_count = -1; 
    cudaError e = cudaGetDeviceCount( &device_count );
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(d_device_number <= device_count, 
            "The number of the device inserted exceeds the set of devices available");

    if( seed_type.c_str()==std::string("same") )   
    	d_seed_type = SEED_TYPE::SAME;
    else if( seed_type.c_str()==std::string("different") )
        d_seed_type = SEED_TYPE::DIFF;
    else if( seed_type.c_str()==std::string("random") ) 
        d_seed_type = SEED_TYPE::RAND;

    // Determine type of tally
    std::string estimator = pl->get<std::string>("estimator",
                                                 "collision");

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

    e = cudaSetDevice( d_device_number );
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

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
    thrust::device_vector<double> rhs(d_N);

    for(int i=0; i<d_N; ++i)
       rhs[i]=b_data[i];

    const double * const rhs_ptr = thrust::raw_pointer_cast(rhs.data());

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

    int block_size;
    int num_blocks;
    int tot_histories;
    
    if(d_use_thread_per_entry)
    {
	//instantiation of as many threads as the number of entries in the solution
    	block_size = std::min(256, d_N);
    	num_blocks = d_N / block_size + 1;
    }
    else 
    {
	//instiantiation of as many threads as the total number of histories
        tot_histories = d_num_histories * d_N;
        block_size = std::min(256,tot_histories);
        num_blocks = tot_histories / block_size + 1;    
    }

    VALIDATE( num_blocks > 0, "The size of the problem is too small" );

    curandState *rng_states;
    cudaError e = cudaMalloc((void **)&rng_states,
        block_size*num_blocks*sizeof(curandState));

    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to allocate memory");

    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    if( d_seed_type==SEED_TYPE::SAME )
    {
        std::cout<<"Same seed instantiated for all the threads"<<std::endl;

    	// Initialize RNG
    	SameSeed seed(d_rng_seed);
	initialize_rng<SameSeed><<<num_blocks,block_size>>>(rng_states,
    	    d_num_curand_calls, seed);
    }
    else if ( d_seed_type==SEED_TYPE::DIFF )
    {
        std::cout<<"Different adjacent seeds instantiated"<<std::endl;

    	DifferentSeed seed( block_size*num_blocks, d_rng_seed);
    	
    	initialize_rng<DifferentSeed><<<num_blocks, block_size>>>(rng_states,  
            d_num_curand_calls, seed);
    }
    else if ( d_seed_type==SEED_TYPE::RAND )
    {
        std::cout<<"Different random seeds instantiated from 0 to "<<
         RAND_MAX<<std::endl;

    	RandomSeed seed( block_size*num_blocks);
 
    	initialize_rng<RandomSeed><<<num_blocks, block_size>>>(rng_states,
            d_num_curand_calls, seed);
    }

    // Check for errors in kernel launch
    e = cudaGetLastError();
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to initialize RNG");
    d_num_curand_calls++;

    if( d_struct==0 )
    {
	if( d_use_thread_per_entry )
        {
            if( d_use_ldg==0 )	
            {   
		    run_forward_monte_carlo2<StandardAccess><<< num_blocks,block_size,sizeof(double)*block_size >>>(d_N,
		    	d_max_history_length, d_weight_cutoff, d_num_histories,
			H,P,W,inds,offsets,coeffs,x_ptr,rhs_ptr, rng_states);  
            }
            else
            {
		    run_forward_monte_carlo2<LDGAccess><<< num_blocks,block_size,sizeof(double)*block_size >>>(d_N,
		    	d_max_history_length, d_weight_cutoff, d_num_histories,
			H,P,W,inds,offsets,coeffs,x_ptr,rhs_ptr, rng_states);  
            }		
        }
	else    
        {
	    int batch_size = 5;
		       
            if( d_use_ldg==0 )		      
            {
		    run_forward_monte_carlo<StandardAccess><<< num_blocks,block_size >>>(d_N,
			d_max_history_length, d_weight_cutoff, d_num_histories, batch_size,
			H,P,W,inds,offsets,coeffs,x_ptr,rhs_ptr, rng_states);            
            }
            else
            {
		    run_forward_monte_carlo<LDGAccess><<< num_blocks,block_size >>>(d_N,
			d_max_history_length, d_weight_cutoff, d_num_histories, batch_size,
			H,P,W,inds,offsets,coeffs,x_ptr,rhs_ptr, rng_states);                   
            }
	}
    }	
    else
    {	
        if( d_use_ldg==0 )
        {
		run_forward_monte_carlo3<StandardAccess><<< num_blocks,block_size >>>( d_N,d_max_history_length,d_weight_cutoff,d_num_histories,data_ptr,
		        coeffs,offsets,x_ptr,rhs_ptr,rng_states  );
        }
        else
        {
		run_forward_monte_carlo3<LDGAccess><<< num_blocks,block_size >>>( d_N,d_max_history_length,d_weight_cutoff,d_num_histories,data_ptr,
		        coeffs,offsets,x_ptr,rhs_ptr,rng_states  );
        }
    }           
    // Check for errors in kernel launch
    e = cudaGetLastError();
    if( cudaSuccess != e )
        std::cout << "Cuda Error: " << cudaGetErrorString(e) << std::endl;

    VALIDATE(cudaSuccess==e,"Failed to execute MC kernel"); 

    // Copy data back to host
    //{
        Teuchos::ArrayRCP<double> x_data = x.getDataNonConst(0);
        thrust::copy(x_vec.begin(),x_vec.end(),x_data.get());
    //}

    // Scale by history count
    for( auto itr= x_data.begin(); itr != x_data.end(); ++itr )
        *itr /= static_cast<double>(d_num_histories); 

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

void ForwardMcCuda::prepareDeviceData(Teuchos::RCP<const MC_Data> mc_data,
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
    		Teuchos::RCP<const MATRIX> H = mc_data->getIterationMatrix();
    		Teuchos::RCP<const MATRIX> P = mc_data->getProbabilityMatrix();
    		Teuchos::RCP<const MATRIX> W = mc_data->getWeightMatrix();

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
            			data_host[count].H = hval_row[j];            
            			data_host[count].P = pval_row[j];
            			data_host[count].W = wval_row[j];
            			data_host[count].inds = ind_row[j];
            			count++;
        		}
        
        		d_offsets[i+1] = d_offsets[i] + ind_row.size();
    		}

    		CHECK( count == d_nnz );

    		thrust::copy( data_host.begin(), data_host.end(), mat_data.begin() );

    		// Copy coefficients into device vector
    		const_scalar_view::HostMirror coeffs_host = Kokkos::create_mirror_view(coeffs);
    		Kokkos::deep_copy(coeffs_host,coeffs);
    		d_coeffs.resize(coeffs.size());
    		thrust::copy(coeffs_host.ptr_on_device(),
                 	coeffs_host.ptr_on_device()+coeffs_host.size(),
                 	d_coeffs.begin());

        }
}


} // namespace alea
