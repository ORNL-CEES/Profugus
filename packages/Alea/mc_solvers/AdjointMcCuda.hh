//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcCuda.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcCuda_hh
#define mc_solver_AdjointMcCuda_hh

#include <string>

#include "MC_Data.hh"

#include "AleaTypedefs.hh"

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include "CudaUtils.hh"
#endif

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class AdjointMcCuda
 * \brief Perform Monte Carlo random walks on linear system with Cuda.
 */
//---------------------------------------------------------------------------//
//

class AdjointMcCuda
{
  public:

    AdjointMcCuda(Teuchos::RCP<const MC_Data>          mc_data,
    		      const const_scalar_view              coeffs,
                  Teuchos::RCP<Teuchos::ParameterList> pl);

    //! Solve problem
    void solve(const MV &b, MV &x);

  private:

#ifdef __CUDACC__
    // Convert from Tpetra matrices to ArrayViews
    void prepareDeviceData(Teuchos::RCP<const MC_Data> mc_data,
            const const_scalar_view coeffs);

    // Build the initial CDF and weights
    void build_initial_distribution(
        Teuchos::ArrayRCP<const double> b,
        thrust::device_vector<double>  &cdf,
        thrust::device_vector<double>  &wt) const;

    // Data for Monte Carlo
    thrust::device_vector<double> d_H;
    thrust::device_vector<double> d_P;
    thrust::device_vector<double> d_W;
    thrust::device_vector<int>    d_inds;
    thrust::device_vector<int>    d_offsets;
    thrust::device_vector<double> d_coeffs;

    thrust::device_vector<device_row_data> mat_data;

    SEED_TYPE d_seed_type;
#endif

    // Vector length
    int d_N;
    int d_nnz;
    int d_num_curand_calls;
    int d_rng_seed;

    // Problem parameters
    int       d_max_history_length;
    bool      d_use_expected_value;
    int       d_num_histories;
    double    d_weight_cutoff;
    bool      d_struct;
    bool      d_use_ldg;
    bool      d_precompute_states;


    enum VERBOSITY {NONE, LOW, HIGH};
    VERBOSITY d_verbosity;
};


} // namespace alea

#endif // mc_solver_AdjointMcCuda_hh

