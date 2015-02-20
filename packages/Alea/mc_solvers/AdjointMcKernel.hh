//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcKernel.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcKernel_hh
#define mc_solver_AdjointMcKernel_hh

#include "MC_Data.hh"
#include "AleaSolver.hh"

#include "Kokkos_View.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class AdjointMcKernel
 * \brief Perform Monte Carlo random walks on linear system.
 *
 * This class performs random walks using the adjoint Monte Carlo algorithm.
 * The interface of this function conforms to the Kokkos "parallel_reduce"
 * functor API to enable automated shared memory parallelism over MC histories.
 */
//---------------------------------------------------------------------------//

class AdjointMcKernel
{
  public:

    // Required typedefs for Kokkos functor API

    //! Type of data kernel will operator on
    typedef SCALAR value_type[];

    // Execution policy and team member types
    typedef Kokkos::TeamPolicy<DEVICE> team_policy;
    typedef team_policy::member_type   team_member;

    // Required public member variable for Kokkos array functor API
    //! Number of entries in a value_type array
    const LO value_count;

    // Convenience typedefs
    typedef Kokkos::View<      SCALAR *, DEVICE> view_type;
    typedef Kokkos::View<const SCALAR *, DEVICE> const_view_type;
    typedef Kokkos::View<const LO     *, DEVICE> const_ord_view;

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;
    typedef typename generator_pool::generator_type generator_type;

    AdjointMcKernel(const const_view_type H,
                    const const_view_type P,
                    const const_view_type W,
                    const const_ord_view  inds,
                    const const_ord_view  offsets,
                    const const_view_type coeffs,
                    const view_type       start_cdf,
                    const view_type       start_wt,
                          int             histories_per_team,
                          bool            use_expected_value,
                          bool            print);

    // Kokkos "parallel_reduce" API functions

    //! \brief Initialize a thread
    KOKKOS_INLINE_FUNCTION
    void init( SCALAR *update ) const;

    //! \brief Compute kernel
    KOKKOS_INLINE_FUNCTION
    void operator()(team_member member, SCALAR *y) const;

    //! \brief Join threads together via a reduce
    KOKKOS_INLINE_FUNCTION
    void join(      volatile SCALAR *update,
              const volatile SCALAR *input) const;

  private:

    KOKKOS_INLINE_FUNCTION
    void getNewRow(const LO        state,
                   const SCALAR * &h_vals,
                   const SCALAR * &p_vals,
                   const SCALAR * &w_vals,
                   const LO     * &inds,
                         LO       &row_length) const;

    KOKKOS_INLINE_FUNCTION
    void tallyContribution(const LO state,
                           const SCALAR wt,
                           const SCALAR * const h_vals,
                           const LO     * const inds,
                           const int            row_length,
                                 SCALAR * const y) const;

    KOKKOS_INLINE_FUNCTION
    LO getNewState(const SCALAR * const  cdf,
                   const LO              cdf_length,
                         generator_type &gen ) const;

    // Data for Monte Carlo
    const const_view_type d_H;
    const const_view_type d_P;
    const const_view_type d_W;
    const const_ord_view  d_inds;
    const const_ord_view  d_offsets;
    const const_view_type d_coeffs;
    const const_view_type d_start_cdf;
    const const_view_type d_start_wt;

    // Kokkos random generator pool
    generator_pool d_rand_pool;

    // Problem parameters
    const int d_histories_per_team;
    const bool d_use_expected_value;
    const bool d_print;
    const int d_max_history_length;

};

} // namespace alea

#include "AdjointMcKernel.i.hh"

#endif // mc_solver_MonteCarloSolver_hh

