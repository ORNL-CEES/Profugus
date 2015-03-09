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

    // lower_bound implementation that can be called from device
    // binary search occurs over range of elements in a Kokkos::View
    template <class view_type>
    KOKKOS_INLINE_FUNCTION
    LO lower_bound(const view_type &v,
                         LO         first,
                         LO         count,
                         SCALAR     val)
    {
        LO current, step;
        while( count > 0 )
        {
            step = count / 2;
            current = first+step;
            if( v(current) < val )
            {
                first = ++current;
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

    InitHistory(scalar_view      randoms,
                const_scalar_view init_cdf,
                const_scalar_view init_wts,
                scalar_view       weights,
                ord_view          states,
                ord_view          starting_inds,
                ord_view          row_lengths,
                ord_view          stages,
                random_ord_view   indices,
                random_ord_view   offsets)
        : d_randoms(randoms)
        , d_init_cdf(init_cdf)
        , d_init_wts(init_wts)
        , d_weights(weights)
        , d_states(states)
        , d_starting_inds(starting_inds)
        , d_row_lengths(row_lengths)
        , d_stages(stages)
        , d_indices(indices)
        , d_offsets(offsets)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        // Perform lower_bound search to get new state
        int state = lower_bound(d_init_cdf,0,d_init_cdf.size(),
                                d_randoms(member));

        d_weights(member)       = d_init_wts(state);
        d_states(member)        = state;
        d_starting_inds(member) = d_offsets(state);
        d_row_lengths(member)   = d_offsets(state+1)-d_offsets(state);
        d_stages(member)        = 0;
    }

  private:

    const scalar_view       d_randoms;
    const const_scalar_view d_init_cdf;
    const const_scalar_view d_init_wts;
    const scalar_view       d_weights;
    const ord_view          d_states;
    const ord_view          d_starting_inds;
    const ord_view          d_row_lengths;
    const ord_view          d_stages;
    const random_ord_view   d_indices;
    const random_ord_view   d_offsets;
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

    StateTransition(scalar_view        randoms,
                    scalar_view        weights,
                    ord_view           states,
                    ord_view           starting_inds,
                    ord_view           row_lengths,
                    ord_view           stages,
                    random_scalar_view P,
                    random_scalar_view W,
                    random_ord_view    indices,
                    random_ord_view    offsets)
        : d_randoms(randoms)
        , d_weights(weights)
        , d_states(states)
        , d_starting_inds(starting_inds)
        , d_row_lengths(row_lengths)
        , d_stages(stages)
        , d_P(P)
        , d_W(W)
        , d_indices(indices)
        , d_offsets(offsets)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        // Perform lower_bound search to get new state
        int new_ind = lower_bound(d_P,d_starting_inds(member),
            d_row_lengths(member),d_randoms(member));

        // Update state
        if( new_ind == d_starting_inds(member)+d_row_lengths(member) )
        {
            d_weights(member) = 0.0;
            d_states(member)  = -1;
        }
        else
        {
            int new_state            = d_indices(new_ind);
            d_weights(member)       *= d_W(new_ind);
            d_states(member)         = new_state;
            d_starting_inds(member)  = d_offsets(new_state);
            d_row_lengths(member)    = d_offsets(new_state+1) -
                                       d_offsets(new_state);
            d_stages(member)++;
        }
    }

  private:

    const scalar_view        d_randoms;
    const scalar_view        d_weights;
    const ord_view           d_states;
    const ord_view           d_starting_inds;
    const ord_view           d_row_lengths;
    const ord_view           d_stages;
    const random_scalar_view d_P;
    const random_scalar_view d_W;
    const random_ord_view    d_indices;
    const random_ord_view    d_offsets;
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

    CollisionTally(scalar_view        weights,
                   ord_view           states,
                   random_scalar_view coeffs,
                   scalar_view        y)
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

    const scalar_view        d_weights;
    const ord_view           d_states;
    const random_scalar_view d_coeffs;
    const scalar_view        d_y;
};

} // end namespace alea

#endif // mc_solvers_MC_Components_hh

//---------------------------------------------------------------------------//
//                 end of MC_Components.hh
//---------------------------------------------------------------------------//
