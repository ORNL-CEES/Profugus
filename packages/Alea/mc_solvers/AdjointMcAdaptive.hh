//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcAdaptive.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcAdaptive_hh
#define mc_solver_AdjointMcAdaptive_hh

#include "MC_Data.hh"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class AdjointMcAdaptive
 * \brief Perform Monte Carlo random walks on linear system.
 *
 * This class performs random walks using the adjoint Monte Carlo algorithm
 * with an adaptive stopping criterion based on variance of the mean.
 */
//---------------------------------------------------------------------------//

class AdjointMcAdaptive
{
  public:

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;
    typedef typename generator_pool::generator_type generator_type;

    AdjointMcAdaptive(Teuchos::RCP<const MC_Data>          mc_data,
                      Teuchos::RCP<Teuchos::ParameterList> pl,
                      generator_pool                       rand_pool);

    //! Solve problem
    void solve(const MV &b, MV &x);

  private:

    // Convert from Tpetra matrices to ArrayViews
    void extractMatrices(Teuchos::RCP<const MC_Data> mc_data);

    // Build the initial CDF and weights
    inline void build_initial_distribution(
        Teuchos::ArrayRCP<const double> b,
        Teuchos::ArrayRCP<double>      &cdf,
        Teuchos::ArrayRCP<double>      &wt) const;

    // Initialize a new history
    inline void initializeHistory(int &state, double &wt,
        Teuchos::ArrayRCP<const double>   start_cdf,
        Teuchos::ArrayRCP<const double>   start_wt,
        Teuchos::ArrayView<const double> &h_row,
        Teuchos::ArrayView<const double> &p_row,
        Teuchos::ArrayView<const double> &w_row,
        Teuchos::ArrayView<const int>    &ind_row);

    // Transition a history to a new state
    inline void getNewState(int &state, double &wt,
        Teuchos::ArrayView<const double> &h_row,
        Teuchos::ArrayView<const double> &p_row,
        Teuchos::ArrayView<const double> &w_row,
        Teuchos::ArrayView<const int>    &ind_row);

    // Add contribution of current history to solution
    inline void tallyContribution(int state, double wt,
        Teuchos::ArrayRCP<double>        x,
        Teuchos::ArrayView<const double> h_row,
        Teuchos::ArrayView<const int>    ind_row) const;

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
    std::string d_verbosity;

    // Kokkos random number generator
    generator_pool d_rand_pool;
    generator_type d_rand_gen;
};

} // namespace alea

#include "AdjointMcAdaptive.i.hh"

#endif // mc_solver_AdjointMcAdaptive_hh

