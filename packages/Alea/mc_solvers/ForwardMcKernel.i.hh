//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ForwardMcKernel.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_ForwardMcKernel_i_hh
#define mc_solver_ForwardMcKernel_i_hh

#include <iterator>
#include <random>

#include "ForwardMcKernel.hh"

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
 * \param x RHS for random walks
 * \param y Solution (output) vector
 * \param rng ThreadedRNG object for generating local random values
 * \param histories_per_state Number of histories to start per state
 * \param print Should debug info be printed?
 */
//---------------------------------------------------------------------------//
ForwardMcKernel::ForwardMcKernel(
        const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > P,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > W,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const LO> > inds,
        const Teuchos::ArrayView<const SCALAR> coeffs,
        const LO local_length,
        const Teuchos::ArrayView<const SCALAR> x,
        const Teuchos::ArrayView<SCALAR> y,
        const int histories_per_state,
        bool print)
  : d_local_length(local_length)
  , d_P(P)
  , d_W(W)
  , d_inds(inds)
  , d_coeffs(coeffs)
  , d_x(x)
  , d_y(y)
  , d_histories_per_state(histories_per_state)
  , d_print(print)
  , d_max_history_length(d_coeffs.size()-1)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void ForwardMcKernel::operator()(member_type member) const
{
    LO     new_ind;
    GO     state;
    Teuchos::ArrayView<const SCALAR> row_cdf, row_wts;
    Teuchos::ArrayView<const GO>     row_inds;

    int myrank = member.team_rank();
    int teamsize = member.team_size();
    int num_my_states = d_local_length / teamsize;
    int num_with_extra = d_local_length % teamsize;
    if( myrank < num_with_extra )
        num_my_states++;

    int mystart;
    if( myrank < num_with_extra )
    {
        mystart = myrank * num_my_states;
    }
    else
    {
        mystart = num_with_extra*(num_my_states+1) +
                  (myrank-num_with_extra)*num_my_states;
    }
    int myend = mystart + num_my_states;

    for( int start_state = mystart; start_state<myend; ++start_state )
    {
        for( int ihist=0; ihist<d_histories_per_state; ++ihist )
        {
            state = start_state;
            SCALAR weight = 1.0;
            SCALAR initial_weight = weight;

            if( d_print )
            {
                std::cout << "Starting history in state " << state
                    << " with initial weight " << initial_weight << std::endl;
            }

            int stage = 0;

            // Transport particle until done
            while(true)
            {
                // Get data and add to tally
                getNewRow(state,row_cdf,row_wts,row_inds);
                d_y[start_state] += d_coeffs[stage]*weight*d_x[state];

                if( stage >= d_max_history_length )
                    break;

                // Get new state index
                new_ind = getNewState(row_cdf,myrank);
                if( new_ind == GO_TRAITS::invalid() )
                    break;

                // Modify weight and update state
                weight *=  row_wts[new_ind];
                state   = row_inds[new_ind];
                stage++;

                if( d_print )
                {
                    std::cout << "Transitioning to state " << state
                        << " with new weight " << weight << std::endl;
                }

            } // while
        } // for ihist
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
void ForwardMcKernel::getNewRow( const GO state,
        Teuchos::ArrayView<const SCALAR> &p_vals,
        Teuchos::ArrayView<const SCALAR> &w_vals,
        Teuchos::ArrayView<const LO>     &inds) const
{
    p_vals = d_P[state];
    w_vals = d_W[state];
    inds   = d_inds[state];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
GO ForwardMcKernel::getNewState(const Teuchos::ArrayView<const SCALAR> cdf,
                                const int thread) const
{
    // Generate random number
    //double rand = d_rng.getRandom(thread);
    double rand = 0.5;

    // Sample cdf to get new state
    Teuchos::ArrayView<const SCALAR>::iterator elem =
        std::lower_bound(cdf.begin(),cdf.end(),rand);

    if( elem == cdf.end() )
        return GO_TRAITS::invalid();

    return elem - cdf.begin();
}

} // namespace alea

#endif // mc_solver_ForwardMcKernel_i_hh
