//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/RNG_Control.cu
 * \author Steven Hamilton
 * \date   Thu Mar 31 09:46:56 2016
 * \brief  RNG_Control class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RNG_Control.cuh"

#include "cuda_utils/Launch_Args.t.cuh"

namespace cuda_mc
{

// Functor to initialize RNG
class RNG_Init
{
  public:

    typedef RNG_Control_DMM::seed_type   seed_type;
    typedef RNG_Control_DMM::RNG_State_t RNG_State_t;

    RNG_Init( RNG_State_t *rngs, seed_type *seeds )
      : d_rngs(rngs)
      , d_seeds(seeds)
    {}

    __device__ void operator()( std::size_t tid ) const
    {
        curand_init(d_seeds[tid],0,0,d_rngs+tid);
    }

  private:

    RNG_State_t *d_rngs;
    seed_type   *d_seeds;
};
        
    

//---------------------------------------------------------------------------//


void RNG_Control_DMM::initialize( int num_streams )
{
    int current_streams = d_rng_states.size();
    int new_streams = num_streams - current_streams;

    // Only initialize new RNG states
    if( new_streams > 0 )
    {
        d_rng_states.resize(num_streams);

        // Create seeds on host
        thrust::host_vector<seed_type> seeds_host(new_streams);
        for( auto &s : seeds_host )
            s = d_dist(d_gen);

        // Copy to device
        thrust::device_vector<seed_type> seeds(seeds_host);

        // Build kernel and launch
        cuda::Launch_Args<cuda::arch::Device> launch_args;
        launch_args.set_num_elements(new_streams);
        RNG_Init init( d_rng_states.data().get() + current_streams,
                       seeds.data().get() );
        cuda::parallel_launch( init, launch_args );

        cudaDeviceSynchronize();
    }
}

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/RNG_Control.cu
//---------------------------------------------------------------------------//
