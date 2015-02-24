//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ForwardMcKernel.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_ForwardMcKernel_hh
#define mc_solver_ForwardMcKernel_hh

#include "MC_Data.hh"
#include "AleaSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Kokkos_Parallel.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class ForwardMcKernel
 * \brief Perform Monte Carlo random walks on linear system.
 *
 * This class performs random walks using the forward Monte Carlo algorithm.
 * The interface of this function conforms to the Kokkos "parallel_reduce"
 * functor API to enable automated shared memory parallelism over MC histories.
 */
//---------------------------------------------------------------------------//

class ForwardMcKernel
{
  public:

    // Required typedefs for Kokkos functor API

    //! Type of device where kernel will be executed
    typedef Kokkos::TeamPolicy<DEVICE> policy_type;
    typedef policy_type::member_type   member_type;

    ForwardMcKernel(const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > P,
                    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > W,
                    const Teuchos::ArrayView<const Teuchos::ArrayView<const LO> > inds,
                    const Teuchos::ArrayView<const SCALAR>  coeffs,
                    const LO local_length,
                    const Teuchos::ArrayView<const SCALAR> x,
                    const Teuchos::ArrayView<SCALAR> y,
                    const int histories_per_state,
                    bool print);

    // Kokkos "parallel_for" API functions

    //! \brief Compute kernel
    KOKKOS_INLINE_FUNCTION
    void operator()(member_type member) const;

  private:

    inline void getNewRow(const GO state,
                          Teuchos::ArrayView<const SCALAR> &p_vals,
                          Teuchos::ArrayView<const SCALAR> &w_vals,
                          Teuchos::ArrayView<const LO>     &inds) const;
    inline GO getNewState(const Teuchos::ArrayView<const SCALAR> cdf,
                          const int thread ) const;

    const LO d_local_length;

    // Warning, these are non reference-counted views
    // The underlying reference-counted objects from which they were
    // created must stay alive for the entire life of this object
    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > d_P;
    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > d_W;
    const Teuchos::ArrayView<const Teuchos::ArrayView<const LO> > d_inds;
    const Teuchos::ArrayView<const SCALAR> d_coeffs;
    const Teuchos::ArrayView<const SCALAR> d_x;
    const Teuchos::ArrayView<SCALAR> d_y;

    //const ThreadedRNG d_rng;

    const int d_histories_per_state;
    const bool d_print;
    const int d_max_history_length;

};

} // namespace alea

#include "ForwardMcKernel.i.hh"

#endif // mc_solver_MonteCarloSolver_hh

