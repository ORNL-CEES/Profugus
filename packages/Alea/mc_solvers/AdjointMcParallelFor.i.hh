//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcParallelFor.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcParallelFor_i_hh
#define mc_solver_AdjointMcParallelFor_i_hh

#include <iterator>
#include <random>
#include <cmath>

#include "AdjointMcParallelFor.hh"
#include "MC_Components.hh"
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
 * \param offsets Starting indices for each matrix row
 * \param coeffs Polynomial coefficients
 * \param pl Problem parameters
 */
//---------------------------------------------------------------------------//
AdjointMcParallelFor::AdjointMcParallelFor(
        const MC_Data_View                  &mc_data,
        const const_scalar_view              coeffs,
        Teuchos::RCP<Teuchos::ParameterList> pl)

  : d_N(mc_data.offsets.size()-1)
  , d_mc_data(mc_data)
  , d_coeffs(coeffs)
  , d_start_cdf("start_cdf",d_N)
  , d_start_wt("start_wt",d_N)
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
void AdjointMcParallelFor::solve(const MV &x, MV &y)
{
    range_policy policy(0,d_num_histories);

    // Build initial probability and weight distributions
    build_initial_distribution(x);

    // Need to get Kokkos view directly, this is silly
    Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);
    const scalar_view y_device("result",d_N);
    const scalar_host_mirror y_mirror =
        Kokkos::create_mirror_view(y_device);

    d_y = y_device;

    // Execute functor
    Kokkos::parallel_for(policy,*this);

    // Copy data back to host
    Kokkos::deep_copy(y_mirror,y_device);
    DEVICE::fence();

    // Apply scale factor
    SCALAR scale_factor = 1.0 / static_cast<SCALAR>(d_num_histories);
    for( LO i=0; i<d_N; ++i )
    {
        y_data[i] = scale_factor*y_mirror(i);
    }

    // Add rhs for expected value
    if( d_use_expected_value )
    {
        scalar_host_mirror coeffs_mirror =
            Kokkos::create_mirror_view(d_coeffs);
        Kokkos::deep_copy(coeffs_mirror,d_coeffs);
        y.update(coeffs_mirror(0),x,1.0);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void AdjointMcParallelFor::operator()(const policy_member &member) const
{
    LO index;
    LO state;
    LO row_length;
    SCALAR weight;

    generator_type rand_gen = d_rand_pool.get_state();

    // Get starting position and weight
    bool valid = initHistory(d_start_cdf,state,index,row_length,weight,rand_gen);
    if( !valid )
    {
        d_rand_pool.free_state(rand_gen);
        return;
    }

    if( std::abs(weight) == 0.0 )
    {
        d_rand_pool.free_state(rand_gen);
        return;
    }
    SCALAR initial_weight = weight;

    if( d_print )
    {
        printf("Thread %i starting history in state %i with initial weight %6.2e\n",
               member,state,initial_weight);
    }

    // Collision estimator starts tallying on zeroth order term
    // Expected value estimator gets this added explicitly at the end
    int stage = 0;
    if( d_use_expected_value )
        stage++;

    // Get data and add to tally
    tallyContribution(state,d_coeffs(stage)*weight);

    // Transport particle until done
    for( ; stage<d_max_history_length; ++stage )
    {
        // Get new state index
        valid = transition(d_mc_data.P,state,index,row_length,weight,rand_gen);
        if( !valid )
        {
            d_rand_pool.free_state(rand_gen);
            return;
        }

        if( d_print )
        {
            printf("Thread %i transitioning to state %i with new weight %6.2e\n",
                   member,state,weight);
        }

        // Get data and add to tally
        tallyContribution(state,d_coeffs(stage)*weight);

    } // while

    d_rand_pool.free_state(rand_gen);
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
void AdjointMcParallelFor::tallyContribution(
        const LO             state,
        const SCALAR         wt ) const
{
    if( d_use_expected_value )
    {
        int start = d_mc_data.offsets(state);
        int row_length = d_mc_data.offsets(state+1)-start;
        if( row_length > 0 )
        {
            for( LO i=0; i<row_length; ++i )
            {
                // P is cdf, we want value of pdf
                //y[inds[i]] += wt*H[i];
                Kokkos::atomic_add(&d_y(d_mc_data.inds[start+i]),
                                   wt*d_mc_data.H[start+i]);
            }
        }
    }
    else
    {
        //y[state] += wt;
        Kokkos::atomic_add(&d_y(state),wt);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a history
 */
//---------------------------------------------------------------------------//
template <class view_type>
bool AdjointMcParallelFor::initHistory(const view_type      &cdf,
                                             LO             &state,
                                             LO             &cdf_start,
                                             LO             &cdf_length,
                                             SCALAR         &weight,
                                             generator_type &gen) const
{
    // Generate random number
    SCALAR rand = Kokkos::rand<generator_type,SCALAR>::draw(gen);

    // Sample cdf to get new state
    // Use local lower_bound implementation, not std library version
    // This allows calling from device
    state = lower_bound(cdf,0,d_N,rand);

    if( state == d_N )
        return false;

    // Get initial state
    weight  =  d_start_wt(state);

    // Get row info
    cdf_start  = d_mc_data.offsets(state);
    cdf_length = d_mc_data.offsets(state+1)-cdf_start;

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
template <class view_type>
bool AdjointMcParallelFor::transition(const view_type      &cdf,
                                            LO             &state,
                                            LO             &cdf_start,
                                            LO             &cdf_length,
                                            SCALAR         &weight,
                                            generator_type &gen) const
{
    // Generate random number
    SCALAR rand = Kokkos::rand<generator_type,SCALAR>::draw(gen);

    // Sample cdf to get new state
    // Use local lower_bound implementation, not std library version
    // This allows calling from device
    LO elem = lower_bound(cdf,cdf_start,cdf_length,rand);

    if( elem == cdf_start+cdf_length )
        return false;

    // Modify weight and update state
    state   =  d_mc_data.inds(elem);
    weight *=  d_mc_data.W(elem);

    // Update current row info
    cdf_start  = d_mc_data.offsets(state);
    cdf_length = d_mc_data.offsets(state+1)-cdf_start;

    return true;
}

//---------------------------------------------------------------------------//
// Build initial cdf and weights
//---------------------------------------------------------------------------//
void AdjointMcParallelFor::build_initial_distribution(const MV &x)
{
    // Build data on host, then explicitly copy to device
    // In future, convert this to a new Kernel to allow building
    //  distributions directly on device if x is allocated there
    scalar_host_mirror start_cdf_host = Kokkos::create_mirror_view(d_start_cdf);
    scalar_host_mirror start_wt_host  = Kokkos::create_mirror_view(d_start_wt);

    Teuchos::ArrayRCP<const SCALAR> x_data = x.getData(0);


    for( LO i=0; i<d_N; ++i )
    {
        start_cdf_host(i) =
            SCALAR_TRAITS::pow(SCALAR_TRAITS::magnitude(x_data[i]),
                               d_start_wt_factor);
    }
    SCALAR pdf_sum = std::accumulate(&start_cdf_host(0),&start_cdf_host(d_N-1)+1,0.0);
    ENSURE( pdf_sum > 0.0 );
    std::transform(&start_cdf_host(0),&start_cdf_host(d_N-1)+1,&start_cdf_host(0),
                   [pdf_sum](SCALAR x){return x/pdf_sum;});
    std::transform(x_data.begin(),x_data.end(),&start_cdf_host(0),
                   &start_wt_host(0),
                   [](SCALAR x, SCALAR y){return y==0.0 ? 0.0 : x/y;});
    std::partial_sum(&start_cdf_host(0),&start_cdf_host(d_N-1)+1,&start_cdf_host(0));
    Kokkos::deep_copy(d_start_cdf,start_cdf_host);
    Kokkos::deep_copy(d_start_wt, start_wt_host);
}


} // namespace alea

#endif // mc_solver_AdjointMcParallelFor_i_hh
