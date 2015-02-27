//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcParallelReduce.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcParallelReduce_i_hh
#define mc_solver_AdjointMcParallelReduce_i_hh

#include <iterator>
#include <random>
#include <cmath>

#include "AdjointMcParallelReduce.hh"
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"

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

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param P Views into entries of probability matrix
 * \param W Views into entries of weight matrix
 * \param inds Views into nonzeros indices
 * \param offsets Starting indices for each matrix row
 * \param coeffs Polynomial coefficients
 * \param pl Problem parameters
 */
//---------------------------------------------------------------------------//
AdjointMcParallelReduce::AdjointMcParallelReduce(const const_view_type                H,
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
  , d_rand_pool(pl->get("random_seed",31891))
  , d_max_history_length(d_coeffs.size()-1)
{
    d_num_histories = pl->get("num_histories",1000);

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
void AdjointMcParallelReduce::solve(const MV &x, MV &y)
{
    // Determine number of histories needed per team
    /*
    int rec_team_size = team_policy::team_size_recommended(*this);
    std::cout << "Recommended team size: " << rec_team_size << std::endl;
    CHECK( rec_team_size > 0 );
    int league_size_req = d_num_histories / rec_team_size;
    std::cout << "Requested league size: " << league_size_req << std::endl;
    CHECK( league_size_req > 0 );
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


    std::cout << "Using team policy with " << num_teams << " teams, "
        << team_size << " threads per team" << std::endl;
    std::cout << "Performing " << total_histories << " total histories, "
        << d_histories_per_team << " per team" << std::endl;
        */

    range_policy policy(0,d_num_histories);

    // Build initial probability and weight distributions
    build_initial_distribution(x);

    // Need to get Kokkos view directly, this is silly
    Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);
    const view_type y_device("result",value_count);
    const view_type::HostMirror y_mirror =
        Kokkos::create_mirror_view(y_device);

    // Execute functor
    Kokkos::parallel_reduce(policy,*this,y_mirror);

    // Apply scale factor
    SCALAR scale_factor = 1.0 / static_cast<SCALAR>(d_num_histories);
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
void AdjointMcParallelReduce::init( SCALAR *update ) const
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
void AdjointMcParallelReduce::operator()(const policy_member &member, SCALAR *y) const
{
    //printf("Executing kernel on team %i of %i, thread %i of %i\n",
    //        member.league_rank(),member.league_size(),member.team_rank(),
    //        member.team_size());
    //printf("Vector address is %p on team %i, thread %i\n",
    //        y,member.league_rank(),member.team_rank());
    LO new_ind;
    LO state;
    const SCALAR * row_h;
    const SCALAR * row_cdf;
    const SCALAR * row_wts;
    const LO     * row_inds;
    int row_length;

    //int team_size = member.team_size();
    //int histories_per_thread = d_histories_per_team / team_size;

    //printf("Getting random number generator on team %i, thread %i\n",
    //        member.league_rank(),member.team_rank());
    generator_type rand_gen = d_rand_pool.get_state();

    int histories_per_thread = 1;
    for( int ihist=0; ihist<histories_per_thread; ++ihist )
    {
        /*
        if( d_print )
        {
            printf("Getting new state on team %i, thread %i\n",
                    member.league_rank(),member.team_rank());
        }
        */

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
            /*
            if( d_print )
            {
                printf("Getting new row on team %i, thread %i\n",
                       member.league_rank(),member.team_rank());
            }
            */
            getNewRow(state,row_h,row_cdf,row_wts,row_inds,row_length);
            /*
            if( d_print )
            {
                printf("Tallying contribution to state %i on team %i, thread %i\n",
                       state,member.league_rank(),member.team_rank());
            }
            */
            tallyContribution(state,d_coeffs(stage)*weight,
                              row_h,row_inds,row_length,y);

            /*
            if( d_print )
            {
                printf("Checking length cutoff on team %i, thread %i\n",
                       member.league_rank(),member.team_rank());
            }
            */
            if( stage >= d_max_history_length )
                break;

            // Get new state index
            /*
            if( d_print )
            {
                printf("Getting new state on team %i, thread %i\n",
                       member.league_rank(),member.team_rank());
            }
            */
            new_ind = getNewState(row_cdf,row_length,rand_gen);
            if( new_ind == -1 )
                break;

            // Modify weight and update state
            weight *=  row_wts[new_ind];
            state   = row_inds[new_ind];
            stage++;

            if( d_print )
            {
                printf("Transitioning to state %i with new weight %6.2e",
                       state,weight);
            }

        } // while
    } // for ihist

    d_rand_pool.free_state(rand_gen);
}

//---------------------------------------------------------------------------//
// Kokkos join
//---------------------------------------------------------------------------//
void AdjointMcParallelReduce::join(      volatile SCALAR *update,
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
void AdjointMcParallelReduce::getNewRow( const LO        state,
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
void AdjointMcParallelReduce::tallyContribution(
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
LO AdjointMcParallelReduce::getNewState(const SCALAR * const  cdf,
                                const LO              cdf_length,
                                      generator_type &gen) const
{
    // Generate random number
    SCALAR rand = Kokkos::rand<generator_type,SCALAR>::draw(gen);

    // Sample cdf to get new state
    // Use local lower_bound implementation, not std library version
    // This allows calling from device
    const SCALAR * const elem = lower_bound(cdf,cdf+cdf_length,rand);

    if( elem == cdf+cdf_length )
        return -1;

    return elem - cdf;
}

//---------------------------------------------------------------------------//
// Build initial cdf and weights
//---------------------------------------------------------------------------//
void AdjointMcParallelReduce::build_initial_distribution(const MV &x)
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

#endif // mc_solver_AdjointMcParallelReduce_i_hh
