//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_solvers/MC_Components.hh
 * \author Steven Hamilton
 * \date   Mon Mar 02 08:35:43 2015
 * \brief  Kernels and utilities for Kokkos Monte Carlo solvers.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_MC_Components_hh
#define mc_solvers_MC_Components_hh

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AleaTypedefs.hh"

namespace alea
{

namespace
{
    // lower_bound implementation that can be called from device
    KOKKOS_INLINE_FUNCTION
    const SCALAR * lower_bound(const SCALAR * first,
                               const SCALAR * last,
                               SCALAR val)
    {
        const SCALAR *it;
        int count, step;
        count = last - first;
        while( count > 0 )
        {
            step = count / 2;
            it = first+step;
            if( *it < val )
            {
                first = ++it;
                count -= step+1;
            }
            else
            {
                count = step;
            }
        }
        return first;
    }
}

//===========================================================================//
/*!
 * \class MC_History
 * \brief Struct for containing data related to MC history for Kokkos kernels
 */
//===========================================================================//

struct MC_History
{
    SCALAR weight;
    LO     state;
    LO     starting_ind;
    LO     row_length;
    LO     stage;
};

//===========================================================================//
/*!
 * \class InitHistory
 * \brief Kernel for initializing histories
 *
 * This class is designed to be used within a Kokkos::parallel_for kernel
 * using a Kokkos::RangePolicy
 */
//===========================================================================//

class InitHistory
{
  public:

    typedef Kokkos::RangePolicy<DEVICE> policy_type;
    typedef policy_type::member_type    policy_member;

    typedef Kokkos::View<      SCALAR     *, DEVICE> view_type;
    typedef Kokkos::View<const SCALAR     *, DEVICE> const_view_type;
    typedef Kokkos::View<      LO         *, DEVICE> ord_view;
    typedef Kokkos::View<const LO         *, DEVICE> const_ord_view;

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;
    typedef typename generator_pool::generator_type generator_type;

    InitHistory(const_view_type init_cdf,
                const_view_type init_wts,
                view_type       weights,
                ord_view        states,
                ord_view        starting_inds,
                ord_view        row_lengths,
                ord_view        stages,
                const_ord_view  indices,
                const_ord_view  offsets)
        : d_init_cdf(init_cdf)
        , d_init_wts(init_wts)
        , d_weights(weights)
        , d_states(states)
        , d_starting_inds(starting_inds)
        , d_row_lengths(row_lengths)
        , d_stages(stages)
        , d_indices(indices)
        , d_offsets(offsets)
        , d_rand_pool(31890)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        // Get random number
        // Is this the best way to use the Kokkos RNG?
        generator_type rand_gen = d_rand_pool.get_state();
        SCALAR rand = Kokkos::rand<generator_type,SCALAR>::draw(rand_gen);
        d_rand_pool.free_state(rand_gen);

        // Perform lower_bound search to get new state
        const SCALAR * const cdf_start = &d_init_cdf(0);
        const SCALAR * const cdf_end   = cdf_start + d_init_cdf.size();
        int state = lower_bound(cdf_start,cdf_end,rand) - cdf_start;

        d_weights(member)       = d_init_wts(state);
        d_states(member)        = state;
        d_starting_inds(member) = d_offsets(state);
        d_row_lengths(member)   = d_offsets(state+1)-d_offsets(state);
        d_stages(member)        = 0;
    }

  private:

    const const_view_type d_init_cdf;
    const const_view_type d_init_wts;
    const view_type       d_weights;
    const ord_view        d_states;
    const ord_view        d_starting_inds;
    const ord_view        d_row_lengths;
    const ord_view        d_stages;
    const const_ord_view  d_indices;
    const const_ord_view  d_offsets;
    generator_pool        d_rand_pool;
};

//===========================================================================//
/*!
 * \class StateTransition
 * \brief Kernel for transitioning a new state
 *
 * This class is designed to be used within a Kokkos::parallel_for kernel
 * using a Kokkos::RangePolicy
 */
//===========================================================================//

class StateTransition
{
  public:

    typedef Kokkos::RangePolicy<DEVICE> policy_type;
    typedef policy_type::member_type    policy_member;

    typedef Kokkos::View<      SCALAR     *, DEVICE> view_type;
    typedef Kokkos::View<const SCALAR     *, DEVICE> const_view_type;
    typedef Kokkos::View<      LO         *, DEVICE> ord_view;
    typedef Kokkos::View<const LO         *, DEVICE> const_ord_view;

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;
    typedef typename generator_pool::generator_type generator_type;

    StateTransition(view_type       weights,
                    ord_view        states,
                    ord_view        starting_inds,
                    ord_view        row_lengths,
                    ord_view        stages,
                    const_view_type P,
                    const_view_type W,
                    const_ord_view  indices,
                    const_ord_view  offsets)
        : d_weights(weights)
        , d_states(states)
        , d_starting_inds(starting_inds)
        , d_row_lengths(row_lengths)
        , d_stages(stages)
        , d_P(P)
        , d_W(W)
        , d_indices(indices)
        , d_offsets(offsets)
        , d_rand_pool(31891)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        // Get random number
        // Is this the best way to use the Kokkos RNG?
        generator_type rand_gen = d_rand_pool.get_state();
        SCALAR rand = Kokkos::rand<generator_type,SCALAR>::draw(rand_gen);
        d_rand_pool.free_state(rand_gen);

        // Perform lower_bound search to get new state
        const SCALAR * const row_start = &d_P(d_starting_inds(member));
        const SCALAR * const row_end   = row_start + d_row_lengths(member);
        int row_index = lower_bound(row_start,row_end,rand) - row_start;

        // Update state
        if( row_index == d_row_lengths(member) )
        {
            d_weights(member) = 0.0;
            d_states(member)  = -1;
        }
        else
        {
            int new_ind   = d_starting_inds(member) + row_index;
            int new_state = d_indices(new_ind);
            d_weights(member)       *= d_W(new_ind);
            d_states(member)         = new_state;
            d_starting_inds(member)  = d_offsets(new_state);
            d_row_lengths(member)    = d_offsets(new_state+1) -
                                       d_offsets(new_state);
            d_stages(member)++;
        }
    }

  private:

    const view_type       d_weights;
    const ord_view        d_states;
    const ord_view        d_starting_inds;
    const ord_view        d_row_lengths;
    const ord_view        d_stages;
    const const_view_type d_P;
    const const_view_type d_W;
    const const_ord_view  d_indices;
    const const_ord_view  d_offsets;
    generator_pool        d_rand_pool;
};

//===========================================================================//
/*!
 * \class CollisionTally
 * \brief Kernel for tallying using collision estimator
 *
 * This class is designed to be used within a Kokkos::parallel_for kernel
 * using a Kokkos::RangePolicy
 */
//===========================================================================//

class CollisionTally
{
  public:

    typedef Kokkos::RangePolicy<DEVICE> policy_type;
    typedef policy_type::member_type    policy_member;

    typedef Kokkos::View<      SCALAR     *, DEVICE> view_type;
    typedef Kokkos::View<const SCALAR     *, DEVICE> const_view_type;
    typedef Kokkos::View<      LO         *, DEVICE> ord_view;
    typedef Kokkos::View<const LO         *, DEVICE> const_ord_view;

    CollisionTally(const_view_type  weights,
                   const_ord_view   states,
                   const_view_type  coeffs,
                   view_type        y)
        : d_weights(weights)
        , d_states(states)
        , d_coeffs(coeffs)
        , d_y(y)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        Kokkos::atomic_add(&d_y(d_states(member)),d_weights(member));
        //d_coeffs(d_histories(member).stage)*d_histories(member).weight);
    }

  private:

    const const_view_type d_weights;
    const const_ord_view  d_states;
    const const_view_type d_coeffs;
    const view_type       d_y;
};

} // end namespace alea

#endif // mc_solvers_MC_Components_hh

//---------------------------------------------------------------------------//
//                 end of MC_Components.hh
//---------------------------------------------------------------------------//
