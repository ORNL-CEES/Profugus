//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcEventKernel.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcEventKernel_hh
#define mc_solver_AdjointMcEventKernel_hh

#include "MC_Data.hh"

#include "Kokkos_View.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class AdjointMcEventKernel
 * \brief Perform Monte Carlo random walks on linear system.
 *
 * This class performs random walks using the adjoint Monte Carlo algorithm.
 * The interface of this function conforms to the Kokkos "parallel_for"
 * functor API to enable automated shared memory parallelism.
 */
//---------------------------------------------------------------------------//

struct NewStateTag {};
struct NewRowTag   {};
struct TallyTag    {};


class AdjointMcEventKernel
{
  public:

    // Required typedefs for Kokkos functor API

    // Execution policy and team member types
    typedef Kokkos::RangePolicy<DEVICE> range_policy;
    typedef range_policy::member_type   policy_member;

    // Convenience typedefs
    typedef Kokkos::View<      SCALAR *,     DEVICE> view_type;
    typedef Kokkos::View<const SCALAR *,     DEVICE> const_view_type;
    typedef Kokkos::View<      LO     *,     DEVICE> ord_view;
    typedef Kokkos::View<const LO     *,     DEVICE> const_ord_view;
    typedef Kokkos::View<      MC_History *, DEVICE> history_view;
    typedef view_type::HostMirror                    host_view_type;

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;
    typedef typename generator_pool::generator_type generator_type;

    AdjointMcEventKernel(const const_view_type                H,
                         const const_view_type                P,
                         const const_view_type                W,
                         const const_ord_view                 inds,
                         const const_ord_view                 offsets,
                         const const_view_type                coeffs,
                         Teuchos::RCP<Teuchos::ParameterList> pl);

    //! Solve problem (this is a host function)
    void solve(const MV &x, MV &y);

  private:

    // Build the initial CDF and weights (host function)
    void build_initial_distribution(const MV &x);

    // Vector length
    int d_N;

    // Data for Monte Carlo
    const const_view_type d_H;
    const const_view_type d_P;
    const const_view_type d_W;
    const const_ord_view  d_inds;
    const const_ord_view  d_offsets;
    const const_view_type d_coeffs;
    const view_type d_start_cdf;
    const view_type d_start_wt;

    view_type d_y;

    // Kokkos random generator pool
    generator_pool d_rand_pool;

    // Problem parameters
    int    d_max_history_length;
    bool   d_use_expected_value;
    bool   d_print;
    int    d_num_histories;
    SCALAR d_start_wt_factor;
};

} // namespace alea

#include "AdjointMcEventKernel.i.hh"

#endif // mc_solver_AdjointMcEventKernel_hh

