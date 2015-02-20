//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcKernel.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcKernel_i_hh
#define mc_solver_AdjointMcKernel_i_hh

#include <iterator>
#include <random>
#include <cmath>

#include "AdjointMcKernel.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param P Views into entries of probability matrix
 * \param W Views into entries of weight matrix
 * \param inds Views into nonzeros indices
 * \param coeffs Polynomial coefficients
 * \param local_length Number of local elements in vector
 * \param start_cdf CDF corresponding to random walk starting locations.
 * \param start_wt  Weights corresponding to random walk starting locations.
 * \param rng ThreadedRNG object for generating local random values
 * \param histories_per_thread Number of histories to be computed
 * \param use_expected_value Should expected value estimator be used?
 * \param print Should debug info be printed?
 */
//---------------------------------------------------------------------------//
AdjointMcKernel::AdjointMcKernel(const const_view_type H,
                                 const const_view_type P,
                                 const const_view_type W,
                                 const const_ord_view  inds,
                                 const const_ord_view  offsets,
                                 const const_view_type coeffs,
                                 const view_type       start_cdf,
                                 const view_type       start_wt,
                                       int             histories_per_team,
                                       bool            use_expected_value,
                                       bool            print)

  : value_count(start_cdf.size())
  , d_H(H)
  , d_P(P)
  , d_W(W)
  , d_inds(inds)
  , d_offsets(offsets)
  , d_coeffs(coeffs)
  , d_start_cdf(start_cdf)
  , d_start_wt(start_wt)
  , d_rand_pool(31891)
  , d_histories_per_team(histories_per_team)
  , d_use_expected_value(use_expected_value)
  , d_print(print)
  , d_max_history_length(d_coeffs.size()-1)
{
}


//---------------------------------------------------------------------------//
// Kokkos init
//---------------------------------------------------------------------------//
void AdjointMcKernel::init( SCALAR *update ) const
{
    for( LO i=0; i<value_count; ++i )
    {
        update[i] = 0.0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void AdjointMcKernel::operator()(team_member member, SCALAR *y) const
{
    LO new_ind;
    LO state;
    const SCALAR * row_h;
    const SCALAR * row_cdf;
    const SCALAR * row_wts;
    const LO     * row_inds;
    int row_length;

    int team_size = member.team_size();
    int histories_per_thread = d_histories_per_team / team_size;

    generator_type rand_gen = d_rand_pool.get_state();

    // What if number of histories per team isn't divisible by team size?
    // Need to adjust history count accordingly

    for( int ihist=0; ihist<histories_per_thread; ++ihist )
    {
        // Get starting position and weight
        state = getNewState(&d_start_cdf(0),value_count,rand_gen);
        if( state == -1 )
            continue;

        if( std::abs(d_start_wt(state)) == 0.0 )
            continue;

        SCALAR weight = d_start_wt(state);
        SCALAR initial_weight = weight;

        if( d_print )
        {
            printf("Starting history in state %i with initial weight %6.2e\n",
                   state,initial_weight);
        }

        // Collision estimator starts tallying on zeroth order term
        // Expected value estimator gets this added explicitly at the end
        int stage = 0;
        if( d_use_expected_value )
            stage++;

        // Transport particle until done
        while(true)
        {
            // Get data and add to tally
            getNewRow(state,row_h,row_cdf,row_wts,row_inds,row_length);
            tallyContribution(state,d_coeffs(stage)*weight,
                              row_h,row_inds,row_length,y);

            if( stage >= d_max_history_length )
                break;

            // Get new state index
            new_ind = getNewState(row_cdf,row_length,rand_gen);
            if( new_ind == -1 )
                break;

            // Modify weight and update state
            weight *=  row_wts[new_ind];
            state   = row_inds[new_ind];
            stage++;

            if( d_print )
            {
                printf("Transitioning to state %i with new weight %6.2e\n",
                       state,weight);
            }

        } // while
    } // for ihist

    d_rand_pool.free_state(rand_gen);
}

//---------------------------------------------------------------------------//
// Kokkos join
//---------------------------------------------------------------------------//
void AdjointMcKernel::join(      volatile SCALAR *update,
                           const volatile SCALAR *input) const
{
    for( LO i=0; i<value_count; ++i )
    {
        update[i] = update[i] + input[i];
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
void AdjointMcKernel::getNewRow( const LO        state,
                                 const SCALAR * &h_vals,
                                 const SCALAR * &p_vals,
                                 const SCALAR * &w_vals,
                                 const LO     * &inds,
                                       LO       &row_length ) const
{
    LO off     = d_offsets(state);
    h_vals     = &d_H(off);
    p_vals     = &d_P(off);
    w_vals     = &d_W(off);
    inds       = &d_inds(off);
    row_length = d_offsets(state+1)-off;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
void AdjointMcKernel::tallyContribution(
        const LO             state,
        const SCALAR         wt,
        const SCALAR * const h_vals,
        const LO     * const inds,
        const int            row_length,
              SCALAR * const y) const
{
    if( d_use_expected_value )
    {
        if( row_length > 0 )
        {
            y[inds[0]] += wt*h_vals[0];
            for( LO i=1; i<row_length; ++i )
            {
                // P is cdf, we want value of pdf
                y[inds[i]] += wt*h_vals[i];
            }
        }
    }
    else
    {
        y[state] += wt;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
LO AdjointMcKernel::getNewState(const SCALAR * const  cdf,
                                const LO              cdf_length,
                                      generator_type &gen) const
{
    // Generate random number
    //SCALAR rand = d_rng.getRandom(thread);
    SCALAR rand = Kokkos::rand<generator_type,SCALAR>::draw(gen);

    // Sample cdf to get new state
    const SCALAR * const elem = std::lower_bound(cdf,cdf+cdf_length,rand);

    if( elem == cdf+cdf_length )
        return -1;

    return elem - cdf;
}

} // namespace alea

#endif // mc_solver_AdjointMcKernel_i_hh
