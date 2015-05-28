//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstMC_Components.cc
 * \author Steven Hamilton
 * \brief  Test of MonteCarloSolver class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <cmath>

#include "../LinearSystem.hh"
#include "../LinearSystemFactory.hh"
#include "../MC_Data.hh"
#include "../MC_Components.hh"
#include "../DeviceTraits.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace alea;

namespace
{
    LO initial_state(LO hist, LO N){return 3*hist % N;}
    SCALAR initial_weight(LO hist){return 1.0+static_cast<SCALAR>(hist);}
    SCALAR my_rand(LO hist)
    {
        return std::fmod(0.05 + 0.1*static_cast<double>(hist),1.0);
    }
}

class TestComponents : public ::testing::TestWithParam<int>
{
  protected:

    void SetUp()
    {
        // Create ParameterList
        d_pl = Teuchos::rcp( new Teuchos::ParameterList() );
        auto mat_pl = Teuchos::sublist(d_pl,"Problem");
        auto mc_pl = Teuchos::sublist(d_pl,"Monte Carlo");
        auto poly_pl = Teuchos::sublist(d_pl,"Polynomial");

        DeviceTraits<DEVICE>::initialize(d_pl);

        d_N = 10;
        mat_pl->set("matrix_type","laplacian");
        mat_pl->set("matrix_size",d_N);

        mc_pl->set("mc_type","adjoint");

        auto system = LinearSystemFactory::buildLinearSystem(d_pl);
        auto A = system->getMatrix();

        Teuchos::RCP<PolynomialBasis> basis(new PolynomialBasis("neumann"));

        Teuchos::RCP<MC_Data> mc_data(new MC_Data(A,basis,d_pl));
        d_mc_data_view = mc_data->createKokkosViews();

        // Get history count from testing framework
        d_num_histories = GetParam();

        // Create history data
        d_hist_data = History_Data(d_num_histories);

        Kokkos::resize(d_randoms,d_num_histories);

        set_initial_state();
    }

    void set_initial_state()
    {
        // Set up initial history data
        auto weight_host = Kokkos::create_mirror_view(d_hist_data.weight);
        auto state_host = Kokkos::create_mirror_view(d_hist_data.state);
        auto starting_ind_host =
            Kokkos::create_mirror_view(d_hist_data.starting_ind);
        auto row_length_host =
            Kokkos::create_mirror_view(d_hist_data.row_length);
        auto offsets_host =
            Kokkos::create_mirror_view(d_mc_data_view.offsets);
        Kokkos::deep_copy(offsets_host,d_mc_data_view.offsets);
        auto rand_host = Kokkos::create_mirror_view(d_randoms);
        for( int hist=0; hist<d_num_histories; ++hist )
        {
            state_host(hist)        = initial_state(hist,d_N);
            weight_host(hist)       = initial_weight(hist);
            starting_ind_host(hist) = offsets_host(state_host(hist));
            row_length_host(hist)   = offsets_host(state_host(hist)+1) -
                                      starting_ind_host(hist);
            rand_host(hist)         = my_rand(hist);
        }
        Kokkos::deep_copy(d_hist_data.weight,weight_host);
        Kokkos::deep_copy(d_hist_data.state,state_host);
        Kokkos::deep_copy(d_hist_data.starting_ind,starting_ind_host);
        Kokkos::deep_copy(d_hist_data.row_length,row_length_host);
        Kokkos::deep_copy(d_randoms,rand_host);
    }

    void TearDown()
    {
        DeviceTraits<DEVICE>::finalize();
    }

    void test_transition()
    {
        // Copy data to host
        auto weight_host = Kokkos::create_mirror_view(d_hist_data.weight);
        Kokkos::deep_copy(weight_host,d_hist_data.weight);
        auto state_host = Kokkos::create_mirror_view(d_hist_data.state);
        Kokkos::deep_copy(state_host,d_hist_data.state);
        auto starting_ind_host =
            Kokkos::create_mirror_view(d_hist_data.starting_ind);
        Kokkos::deep_copy(starting_ind_host,d_hist_data.starting_ind);
        auto row_length_host =
            Kokkos::create_mirror_view(d_hist_data.row_length);
        Kokkos::deep_copy(row_length_host,d_hist_data.row_length);
        auto offsets_host =
            Kokkos::create_mirror_view(d_mc_data_view.offsets);
        Kokkos::deep_copy(offsets_host,d_mc_data_view.offsets);
        auto rand_host = Kokkos::create_mirror_view(d_randoms);
        Kokkos::deep_copy(rand_host,d_randoms);

        // Currently we check that each individual history undergoes the
        // correct transition, in the future it may be beneficial to allow
        // transition kernels to reorder the histories.  To allow that we
        // would need to refactor this test to make sure that any history
        // is in the expected state rather than requiring it of a particular
        // history.
        for( int hist=0; hist<d_num_histories; ++hist )
        {
            // Determine original state of history
            LO     init_state = initial_state(hist,d_N);
            SCALAR init_wt    = initial_weight(hist);
            SCALAR this_rand = my_rand(hist);

            // Check if transition was correct
            // States on edges move away from boundary and have weight reduced
            if( init_state == 0 )
            {
                EXPECT_EQ(state_host(hist),1);
                EXPECT_DOUBLE_EQ(init_wt/2.0,weight_host(hist));
            }
            else if( init_state == d_N-1 )
            {
                EXPECT_EQ(state_host(hist),d_N-2);
                EXPECT_DOUBLE_EQ(init_wt/2.0,weight_host(hist));
            }
            // States just off boundary have different probability distribution
            else if( init_state == 1 )
            {
                if( this_rand < 0.4 )
                {
                    EXPECT_EQ(state_host(hist),0);
                }
                else
                {
                    EXPECT_EQ(state_host(hist),2);
                }
                EXPECT_DOUBLE_EQ(init_wt*5.0/6.0,weight_host(hist));
            }
            else if( init_state == d_N-2 )
            {
                if( this_rand < 0.6 )
                {
                    EXPECT_EQ(state_host(hist),d_N-3);
                }
                else
                {
                    EXPECT_EQ(state_host(hist),d_N-1);
                }
                EXPECT_DOUBLE_EQ(init_wt*5.0/6.0,weight_host(hist));
            }
            // States in interior move left or right and have weight unchanged
            else
            {
                if( this_rand < 0.5 )
                    EXPECT_EQ(init_state-1,state_host(hist));
                else
                    EXPECT_EQ(init_state+1,state_host(hist));

                EXPECT_DOUBLE_EQ(init_wt,weight_host(hist));
            }

            // Not all kernels use the starting index and row length so we
            // can't verify them
        }
    }

    Teuchos::RCP<Teuchos::ParameterList> d_pl;
    MC_Data_View d_mc_data_view;
    History_Data d_hist_data;
    scalar_view d_randoms;
    int d_N;
    int d_num_histories;
};

// Run each test with several history counts
INSTANTIATE_TEST_CASE_P(Default,TestComponents,
                        ::testing::Values(10,100,1000,10000));

TEST_P(TestComponents,StandardTransition)
{
    // Build policy and kernel and execute
    Kokkos::RangePolicy<DEVICE> policy(0,d_randoms.size());
    StateTransition transition(d_randoms,d_hist_data,d_mc_data_view);
    Kokkos::parallel_for(policy,transition);
    DEVICE::fence();

    // Test correctness of result
    test_transition();
}

TEST_P(TestComponents,BinnedTransition)
{
    ord_view old_states("old_states",d_num_histories);

    // Build policy and kernel and execute
    BinnedStateTransition transition(
        d_randoms,old_states,d_hist_data,d_mc_data_view,d_N);

    int team_size = Kokkos::TeamPolicy<DEVICE>::team_size_max(transition);

    // Execute kernel for different league sizes
    std::vector<int> league_size = {1,4,16,64};
    for( int i=0; i<league_size.size(); ++i )
    {
        set_initial_state();

        Kokkos::deep_copy(old_states,d_hist_data.state);
        Kokkos::TeamPolicy<DEVICE> policy(league_size[i],team_size);
        Kokkos::parallel_for(policy,transition);
        DEVICE::fence();

        test_transition();
    }
}

TEST_P(TestComponents,SharedMemTransition)
{
    ord_view old_states("old_states",d_num_histories);

    std::vector<int> num_shared_values = {4, 8, 16};
    for( int i=0; i<num_shared_values.size(); ++i )
    {
        // Build policy and kernel and execute
        SharedMemTransition transition(d_randoms,old_states,d_hist_data,
            d_mc_data_view,d_N,num_shared_values[i]);

        int team_size = Kokkos::TeamPolicy<DEVICE>::team_size_max(transition);

        // Execute kernel for different league sizes
        int league_size = transition.num_blocks();
        Kokkos::TeamPolicy<DEVICE> policy(league_size,team_size);

        set_initial_state();

        Kokkos::deep_copy(old_states,d_hist_data.state);
        Kokkos::parallel_for(policy,transition);
        DEVICE::fence();

        test_transition();
    }
}
