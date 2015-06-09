//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ForwardMcAdaptive.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_ForwardMcAdaptive_hh
#define mc_solver_ForwardMcAdaptive_hh

#include "MC_Data.hh"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class ForwardMcAdaptive
 * \brief Perform Monte Carlo random walks on linear system.
 *
 * This class performs random walks using the forward Monte Carlo algorithm
 * with an adaptive stopping criterion based on variance of the mean.
 */
//---------------------------------------------------------------------------//

class ForwardMcAdaptive
{
  public:

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;
    typedef typename generator_pool::generator_type generator_type;

    ForwardMcAdaptive(Teuchos::RCP<const MC_Data>          mc_data,
    		      const const_scalar_view              coeffs,   // modified by Max
                      Teuchos::RCP<Teuchos::ParameterList> pl,
                      generator_pool                       rand_pool);

    //! Solve problem
    void solve(const MV &b, MV &x);

  private:

    // Convert from Tpetra matrices to ArrayViews
    void extractMatrices(Teuchos::RCP<const MC_Data> mc_data);

    // Transition a history to a new state
    inline void getNewState(int &state, double &wt,
        Teuchos::ArrayView<const double> &h_row,
        Teuchos::ArrayView<const double> &p_row,
        Teuchos::ArrayView<const double> &w_row,
        Teuchos::ArrayView<const int>    &ind_row);

    // Add contribution of current history to solution
    inline void tallyContribution(double wt,double& x);

    // Data for Monte Carlo
    Teuchos::ArrayRCP<Teuchos::ArrayView<const double> > d_H;
    Teuchos::ArrayRCP<Teuchos::ArrayView<const double> > d_P;
    Teuchos::ArrayRCP<Teuchos::ArrayView<const double> > d_W;
    Teuchos::ArrayRCP<Teuchos::ArrayView<const int> >    d_ind;

    // Vector length
    int d_N;

    // Problem parameters
    int    d_max_history_length;
    bool   d_use_expected_value;
    int    d_max_num_histories;
    double d_weight_cutoff;
    int    d_batch_size;
    double d_tolerance;
    const const_scalar_view d_coeffs; // Modified by Max

    enum VERBOSITY {NONE, LOW, HIGH};
    VERBOSITY d_verbosity;

    // Kokkos random number generator
    generator_pool d_rand_pool;
    generator_type d_rand_gen;
};

} // namespace alea

#include "ForwardMcAdaptive.i.hh"

#endif // mc_solver_ForwardMcAdaptive_hh

