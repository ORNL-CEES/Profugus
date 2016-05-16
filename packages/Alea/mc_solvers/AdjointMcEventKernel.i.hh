//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/AdjointMcEventKernel.i.hh
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_AdjointMcEventKernel_i_hh
#define Alea_mc_solvers_AdjointMcEventKernel_i_hh

#include <iterator>
#include <random>
#include <cmath>

#if defined(__CUDACC__)
#include <thrust/device_ptr.h>
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>
#endif

#include "AdjointMcEventKernel.hh"
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
AdjointMcEventKernel::AdjointMcEventKernel(
        const MC_Data_View                  &mc_data,
        const const_scalar_view              coeffs,
        Teuchos::RCP<Teuchos::ParameterList> pl,
        generator_pool                       pool)
  : d_N(mc_data.offsets.size()-1)
  , d_mc_data(mc_data)
  , d_coeffs(coeffs)
  , d_start_cdf("start_cdf",d_N)
  , d_start_wt("start_wt",d_N)
  , d_rand_pool(pool)
  , d_max_history_length(d_coeffs.size()-1)
  , d_pl(pl)
{
    d_num_histories = pl->get("num_histories",1000);

    d_num_batches = pl->get("num_batches",1);
    d_histories_batch = d_num_histories / d_num_batches;

    if( d_num_batches > 1 )
    {
        std::cout << "Performing " << d_num_histories << " per iteration in "
            << d_num_batches << " batches with " << d_histories_batch
            << " histories per batch" << std::endl;
    }

    // Determine type of tally
    std::string estimator = pl->get<std::string>("estimator",
                                                 "expected_value");
    VALIDATE(estimator == "collision" ||
             estimator == "expected_value",
             "Only collision and expected_value estimators are available.");
    d_use_expected_value = (estimator == "expected_value");
    VALIDATE(!d_use_expected_value,
        "Expected value not available in event kernel yet.");

    // Power factor for initial probability distribution
    d_start_wt_factor = pl->get<SCALAR>("start_weight_factor",1.0);

    // Should we print anything to screen
    std::string verb = profugus::lower(pl->get("verbosity","low"));
    d_print = (verb == "high");

    std::string transition =
        profugus::lower(pl->get("transition_type","standard"));
    if( transition == "standard" )
        d_transition_type = STANDARD;
    else if( transition == "binned" )
        d_transition_type = BINNED;
    else if( transition == "shared_mem" )
        d_transition_type = SHARED_MEM;
    else
    {
        std::stringstream ss;
        ss << "Unknown transition kernel type: " << transition << std::endl;
        VALIDATE(false,ss.str());
    }
}

//---------------------------------------------------------------------------//
// Solve problem using Monte Carlo
//---------------------------------------------------------------------------//
void AdjointMcEventKernel::solve(const MV &x, MV &y)
{
    // Build initial probability and weight distributions
    build_initial_distribution(x);

    // Need to get Kokkos view directly, this is silly
    const scalar_view  y_device( "result", d_N);

    // Zero out y
    ZeroVector      zero_y(y_device);
    Kokkos::parallel_for(d_N,zero_y);

    solve_impl(y_device);

    // Get view of data on host
    const scalar_host_mirror y_mirror =
        Kokkos::create_mirror_view(y_device);
    SCALAR scale_factor = 1.0 / static_cast<SCALAR>(d_num_histories);

    // Copy data back to host
    Kokkos::deep_copy(y_mirror,y_device);

    // Apply scale factor
    // Here we force a copy of the data in y to the CPU because we can't
    // currently access it directly on the GPU due to Tpetra limitations.
    // The surrounding braces force a copy of the data back to the GPU before
    // moving to the expected value update.
    {
        Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);
        for( LO i=0; i<d_N; ++i )
        {
            y_data[i] = scale_factor*y_mirror(i);
        }
    }

    // Add rhs for expected value
    if( d_use_expected_value )
    {
        y.update(d_coeffs(0),x,1.0);
    }
}

//---------------------------------------------------------------------------//
// Sort history data by state
//---------------------------------------------------------------------------//
template <class device>
void AdjointMcEventKernel::sort_by_state(const History_Data &data) const
{
    const ord_view state = data.state;
    const scalar_view weight = data.weight;
    const ord_view starting_ind = data.starting_ind;
    const ord_view row_length = data.row_length;
    int N = state.size();

    // Pack data into tuple
    std::vector<std::tuple<LO,SCALAR,LO,LO> > v(N);
    for( LO i=0; i<N; ++i )
        v[i] = std::make_tuple(
            state(i),weight(i),starting_ind(i),row_length(i));

    // Sort
    std::sort(v.begin(),v.end());

    // Put back into Views
    for( LO i=0; i<N; ++i )
    {
        state(i)        = std::get<0>(v[i]);
        weight(i)       = std::get<1>(v[i]);
        starting_ind(i) = std::get<2>(v[i]);
        row_length(i)   = std::get<3>(v[i]);
    }
}

#ifdef KOKKOS_HAVE_CUDA
//---------------------------------------------------------------------------//
// Sort history data by state -- Cuda specialization using thrust
//---------------------------------------------------------------------------//
template <>
void AdjointMcEventKernel::sort_by_state<Kokkos::Cuda>(
    const History_Data &data) const
{
    // Get Kokkos Views
    const ord_view state = data.state;
    const scalar_view weight = data.weight;
    const ord_view starting_ind = data.starting_ind;
    const ord_view row_length = data.row_length;
    int N = state.size();

    // Convert raw pointers to thrust device pointers
    auto state_ptr = thrust::device_pointer_cast(state.ptr_on_device());
    auto weight_ptr = thrust::device_pointer_cast(weight.ptr_on_device());
    auto starting_ind_ptr = thrust::device_pointer_cast(
        starting_ind.ptr_on_device());
    auto row_length_ptr = thrust::device_pointer_cast(
        row_length.ptr_on_device());

    // Sort using zip iterator
    auto tuple_begin = thrust::make_zip_iterator( thrust::make_tuple(
        state_ptr,weight_ptr,starting_ind_ptr,row_length_ptr));
    auto tuple_end = thrust::make_zip_iterator( thrust::make_tuple(
        state_ptr+N,weight_ptr+N,starting_ind_ptr+N,row_length_ptr+N));
    thrust::sort(tuple_begin,tuple_end);
}
#endif

//---------------------------------------------------------------------------//
// Solve implementation for global memory kernel
//---------------------------------------------------------------------------//
void AdjointMcEventKernel::solve_impl(const scalar_view &y_device) const
{
    const scalar_view  randoms(  "randoms",d_histories_batch);
    const History_Data hist_data(d_histories_batch);

    // The binning process for the binned and shared memory kernels need
    // access to the old list of states to avoid a race condition
    ord_view old_states("old_states");
    if( d_transition_type==BINNED ||
        d_transition_type==SHARED_MEM )
    {
        Kokkos::resize(old_states,d_histories_batch);
    }

    // Build kernels
    InitHistory     init_history(randoms,d_start_cdf,d_start_wt,hist_data,
                                 d_mc_data);
    StateTransition transition(randoms,hist_data,d_mc_data);
    CollisionTally  coll_tally(hist_data,d_coeffs,y_device);
    BinnedStateTransition binned_transition(randoms,old_states,hist_data,
        d_mc_data,y_device.size());

    int num_shared_values = d_pl->get("num_shared_values",512);
    SharedMemTransition shared_mem_transition(randoms,old_states,hist_data,
        d_mc_data,y_device.size(),num_shared_values);

    // Create policy
    range_policy policy(0,d_histories_batch);

    // Determine league size -- shared memory kernel requires a particular
    // league size, the other kernels can use any size determined by an
    // input parameter
    int league_size;
    if( d_transition_type == SHARED_MEM )
        league_size = shared_mem_transition.num_blocks();
    else
        league_size = d_pl->get<int>("league_size",128);

    int team_size = d_pl->get<int>("team_size",
        team_policy::team_size_max(binned_transition));
    team_policy team_pol(league_size,team_size);

    std::cout << "Team policy has league size of " << team_pol.league_size()
        << " and team size of " << team_pol.team_size() << std::endl;

    for( int batch=0; batch<d_num_batches; ++batch )
    {
        // Get initial state and tally
        Kokkos::fill_random(randoms,d_rand_pool,1.0);
        Kokkos::parallel_for(policy,init_history);
        Kokkos::parallel_for(policy,coll_tally);

        // Loop over history length (start at 1)
        for( int step=1; step<=d_max_history_length; ++step )
        {
            //sort_by_state<DEVICE>(hist_data);
            Kokkos::fill_random(randoms,d_rand_pool,1.0);
            switch( d_transition_type )
            {
                case STANDARD:
                    Kokkos::parallel_for(policy,transition);
                    break;
                case BINNED:
                    Kokkos::deep_copy(old_states,hist_data.state);
                    Kokkos::parallel_for(team_pol,binned_transition);
                    break;
                case SHARED_MEM:
                    Kokkos::deep_copy(old_states,hist_data.state);
                    Kokkos::parallel_for(team_pol,shared_mem_transition);
                    break;
            };
            Kokkos::parallel_for(policy,coll_tally);
        }
    }
}

//---------------------------------------------------------------------------//
// Build initial cdf and weights
//---------------------------------------------------------------------------//
void AdjointMcEventKernel::build_initial_distribution(const MV &x)
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

#endif // Alea_mc_solvers_AdjointMcEventKernel_i_hh
