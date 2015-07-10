//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils.hh
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_CudaUtils_hh
#define mc_solver_CudaUtils_hh

namespace alea
{

// lower_bound implementation that can be called from device
__device__ const double * lower_bound(const double * first,
        const double * last,
        double   val);
        
__device__ double atomicAdd(double* address, double val);        
        
__device__ void getNewState(int &state, double &wt,
        const double * const P,
        const double * const W,
        const int    * const inds,
        const int    * const offsets,
              curandState   *rng_state );
              
__device__ void getNewState2(int &state, double &wt,
        const double * const P,
        const double * const W,
        const int    * const inds,
        const int    * const offsets,
              double   &rand );                      


__global__ void initialize_rng(curandState *state, int seed, int offset);

__global__ void initialize_rng2(curandState *state, int*seed, int offset)

        
}


#endif // mc_solver_CudaUtils_hh
