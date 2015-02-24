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
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"

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
AdjointMcKernel::AdjointMcKernel(const const_view_type                H,
                                 const const_view_type                P,
                                 const const_view_type                W,
                                 const const_ord_view                 inds,
                                 const const_ord_view                 offsets,
                                 const const_view_type                coeffs,
                                 Teuchos::RCP<Teuchos::ParameterList> pl)

  : value_count(offsets.size()-1)
  , d_H(H)
  , d_P(P)
  , d_W(W)
  , d_inds(inds)
  , d_offsets(offsets)
  , d_coeffs(coeffs)
  , d_start_cdf("start_cdf",value_count)
  , d_start_wt("start_wt",value_count)
  , d_max_history_length(d_coeffs.size()-1)
{
    d_num_histories = pl->get("num_histories",1000);

    // Set up RNG
    int rand_seed = pl->get("random_seed",31891);
    d_rand_pool.init(rand_seed,DEVICE::max_hardware_threads());

    // Determine type of tally
    std::string estimator = pl->get<std::string>("estimator",
                                                 "expected_value");
    VALIDATE(estimator == "collision" ||
             estimator == "expected_value",
             "Only collision and expected_value estimators are available.");
    d_use_expected_value = (estimator == "expected_value");

    // Power factor for initial probability distribution
    d_start_wt_factor = pl->get<SCALAR>("start_weight_factor",1.0);

    // Should we print anything to screen
    std::string verb = profugus::to_lower(pl->get("verbosity","low"));
    d_print = (verb == "high");
}

//---------------------------------------------------------------------------//
// Solve problem using Monte Carlo
//---------------------------------------------------------------------------//
void AdjointMcKernel::solve(const MV &x, MV &y)
{
    // Determine number of histories needed per team
    int rec_team_size = team_policy::team_size_recommended(*this);
    int league_size_req = d_num_histories / rec_team_size;
    team_policy policy(league_size_req,rec_team_size);
    int num_teams = policy.league_size();
    int team_size = policy.team_size();
    int num_threads = num_teams * team_size;
    int total_histories = d_num_histories;
    if( d_num_histories % num_threads != 0 )
    {
        total_histories = (total_histories/num_threads + 1)*num_threads;
        ADD_WARNING("Requested number of histories (" << d_num_histories
           << ") is not divisible by number of threads ("
           << num_threads << "), number of histories is being increased to "
           << total_histories << std::endl);
    }
    d_histories_per_team = total_histories / num_teams;

    // Build initial probability and weight distributions
    build_initial_distribution(x);

    // Need to get Kokkos view directly, this is silly
    Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);
    const view_type y_device("result",value_count);
    const view_type::HostMirror y_mirror =
        Kokkos::create_mirror_view(y_device);

    // Execute functor
    Kokkos::parallel_reduce(policy,*this,y_mirror);

    // Copy data back to host, this shouldn't need to happen
    Kokkos::deep_copy(y_mirror,y_device);

    // Apply scale factor
    SCALAR scale_factor = 1.0 / static_cast<SCALAR>(total_histories);
    for( LO i=0; i<value_count; ++i )
    {
        y_data[i] = scale_factor*y_mirror(i);
    }

    // Add rhs for expected value
    if( d_use_expected_value )
    {
        y.update(d_coeffs(0),x,1.0);
    }
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

// Build initial cdf and weights
void AdjointMcKernel::build_initial_distribution(const MV &x)
{
    // Build data on host, then explicitly copy to device
    // In future, convert this to a new Kernel to allow building
    //  distributions directly on device if x is allocated there
    host_view_type start_cdf_host = Kokkos::create_mirror_view(d_start_cdf);
    host_view_type start_wt_host  = Kokkos::create_mirror_view(d_start_wt);

    Teuchos::ArrayRCP<const SCALAR> x_data = x.getData(0);

    int N = value_count;

    for( LO i=0; i<N; ++i )
    {
        start_cdf_host(i) =
            SCALAR_TRAITS::pow(SCALAR_TRAITS::magnitude(x_data[i]),
                               d_start_wt_factor);
    }
    SCALAR pdf_sum = std::accumulate(&start_cdf_host(0),&start_cdf_host(N-1)+1,0.0);
    ENSURE( pdf_sum > 0.0 );
    std::transform(&start_cdf_host(0),&start_cdf_host(N-1)+1,&start_cdf_host(0),
                   [pdf_sum](SCALAR x){return x/pdf_sum;});
    std::transform(x_data.begin(),x_data.end(),&start_cdf_host(0),
                   &start_wt_host(0),
                   [](SCALAR x, SCALAR y){return y==0.0 ? 0.0 : x/y;});
    std::partial_sum(&start_cdf_host(0),&start_cdf_host(N-1)+1,&start_cdf_host(0));
    Kokkos::deep_copy(d_start_cdf,start_cdf_host);
    Kokkos::deep_copy(d_start_wt, start_wt_host);
}


} // namespace alea

#endif // mc_solver_AdjointMcKernel_i_hh
