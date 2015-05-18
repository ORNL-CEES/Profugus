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
 * \class History_Data
 * \brief Struct for containing data related to MC history for Kokkos kernels
 */
//===========================================================================//
struct History_Data
{
    History_Data(){}
    History_Data(int N)
        : weight("weight",N)
        , state("state",N)
        , starting_ind("starting_ind",N)
        , row_length("row_length",N)
        , stage("stage",N)
    {}

    scalar_view weight;
    ord_view    state;
    ord_view    starting_ind;
    ord_view    row_length;
    ord_view    stage;
};

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
 * \class ZeroVector
 * \brief Zeros out a vector
 *
 * Doesn't Kokkos have something built-in for this?
 */
//===========================================================================//

class ZeroVector
{
  public:

    typedef Kokkos::RangePolicy<DEVICE> policy_type;
    typedef policy_type::member_type    policy_member;

    ZeroVector(scalar_view v)
        : d_v(v)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const
    {
        d_v(i) = 0.0;
    }

  private:

    const scalar_view d_v;
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

    InitHistory(scalar_view         randoms,
                const_scalar_view   init_cdf,
                const_scalar_view   init_wts,
                const History_Data &hist_data,
                const MC_Data_View &mc_data)
        : d_randoms(randoms)
        , d_init_cdf(init_cdf)
        , d_init_wts(init_wts)
        , d_hist_data(hist_data)
        , d_mc_data(mc_data)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        // Perform lower_bound search to get new state
        int state = lower_bound(d_init_cdf,0,d_init_cdf.size(),
                                d_randoms(member));

        d_hist_data.weight(member)       = d_init_wts(state);
        d_hist_data.state(member)        = state;
        d_hist_data.starting_ind(member) = d_mc_data.offsets(state);
        d_hist_data.row_length(member)   = d_mc_data.offsets(state+1) -
                                           d_mc_data.offsets(state);
        d_hist_data.stage(member)        = 0;
    }

  private:

    const random_scalar_view   d_randoms;
    const random_scalar_view   d_init_cdf;
    const random_scalar_view   d_init_wts;
    const History_Data         d_hist_data;
    const MC_Data_Texture_View d_mc_data;
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

    StateTransition(scalar_view         randoms,
                    const History_Data &hist_data,
                    const MC_Data_View &mc_data)
        : d_randoms(randoms)
        , d_hist_data(hist_data)
        , d_mc_data(mc_data)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        // Perform lower_bound search to get new state
        int new_ind = lower_bound(d_mc_data.P,d_hist_data.starting_ind(member),
            d_hist_data.row_length(member),d_randoms(member));

        // Update state
        if( new_ind == d_hist_data.starting_ind(member) +
                       d_hist_data.row_length(member) )
        {
            d_hist_data.weight(member) = 0.0;
            d_hist_data.state(member)  = -1;
        }
        else
        {
            int new_state                    = d_mc_data.inds(new_ind);
            d_hist_data.weight(member)      *= d_mc_data.W(new_ind);
            d_hist_data.state(member)        = new_state;
            d_hist_data.starting_ind(member) = d_mc_data.offsets(new_state);
            d_hist_data.row_length(member)   = d_mc_data.offsets(new_state+1) -
                                               d_mc_data.offsets(new_state);
            d_hist_data.stage(member)++;
        }
    }

  private:

    const random_scalar_view   d_randoms;
    const History_Data         d_hist_data;
    const MC_Data_Texture_View d_mc_data;
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

    CollisionTally(const History_Data &hist_data,
                   random_scalar_view coeffs,
                   scalar_view        y)
        : d_hist_data(hist_data)
        , d_coeffs(coeffs)
        , d_y(y)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(policy_member member) const
    {
        Kokkos::atomic_add(&d_y(d_hist_data.state(member)),d_hist_data.weight(member));
        //d_coeffs(d_histories(member).stage)*d_histories(member).weight);
    }

  private:

    const History_Data       d_hist_data;
    const random_scalar_view d_coeffs;
    const scalar_view        d_y;
};

//===========================================================================//
/*!
 * \class TeamEventKernel
 * \brief Team-based kernel for event MC
 *
 * This class is designed to be used within a Kokkos::parallel_for kernel
 * using a Kokkos::TeamPolicy
 */
//===========================================================================//

class TeamEventKernel
{
  public:

    typedef Kokkos::TeamPolicy<DEVICE> policy_type;
    typedef policy_type::member_type   policy_member;

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>      generator_pool;

    unsigned team_shmem_size(int team_size) const
    {
        //return team_size*(sizeof(SCALAR)+4*sizeof(LO));
        int histories_team = d_num_histories / d_league_size;
        return histories_team * (sizeof(SCALAR)+3*sizeof(LO));
    }

    TeamEventKernel( int max_history_length,
                     int num_histories,
                     const_scalar_view   start_cdf,
                     const_scalar_view   start_wt,
                     const MC_Data_View &mc_data,
                     scalar_view         y,
                     random_scalar_view  coeffs )
        : d_max_history_length(max_history_length)
        , d_num_histories(num_histories)
        , d_start_cdf(start_cdf)
        , d_start_wt(start_wt)
        , d_mc_data(mc_data)
        , d_y(y)
        , d_coeffs(coeffs)
    {
    }

    void set_league_size(int league_size){d_league_size = league_size;}

    void set_randoms(scalar_view_2d randoms){d_randoms = randoms;}

    KOKKOS_INLINE_FUNCTION
    void operator()(const policy_member &member) const
    {
        int team_size = member.team_size();
        int histories_team = d_num_histories / member.league_size();
        int team_offset = histories_team*member.league_rank();

        shared_scalar_view weight(       member.team_shmem(), histories_team);
        shared_ord_view    state(        member.team_shmem(), histories_team);
        shared_ord_view    starting_ind( member.team_shmem(), histories_team);
        shared_ord_view    row_length(   member.team_shmem(), histories_team);

        // Bring member variables into local scope to allow lambda capture
        int stage = 0;

        // Init kernel
        auto init = [=](const int &i)
        {
            // Perform lower_bound search to get new state
            int init_state = lower_bound(d_start_cdf,0,d_start_cdf.size(),
                d_randoms(i+team_offset,stage));

            weight(i)       = d_start_wt(init_state);
            state(i)        = init_state;
            starting_ind(i) = d_mc_data.offsets(init_state);
            row_length(i)   = d_mc_data.offsets(init_state+1) -
                              d_mc_data.offsets(init_state);
        };

        auto transition = [&](const int &i)
        {
            // Perform lower_bound search to get new state
            int new_ind = lower_bound(d_mc_data.P,starting_ind(i),
                row_length(i),d_randoms(i+team_offset,stage));

            // Update state
            if( new_ind == starting_ind(i) + row_length(i) )
            {
                weight(i) = 0.0;
                state(i)  = -1;
            }
            else
            {
                int new_state   = d_mc_data.inds(new_ind);
                weight(i)      *= d_mc_data.W(new_ind);
                state(i)        = new_state;
                starting_ind(i) = d_mc_data.offsets(new_state);
                row_length(i)   = d_mc_data.offsets(new_state+1) -
                                  d_mc_data.offsets(new_state);
            }
        };

        // Tally kernel
        auto tally = [=](const int &i)
        {
            Kokkos::atomic_add(&d_y(state(i)),weight(i));
        };

        // Do history initialization and initial tally
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(member,histories_team),init);
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(member,histories_team),tally);

        // Process remaining steps for all histories
        for( int step=0; step<d_max_history_length; ++step )
        {
            stage++;
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(member,histories_team),transition);
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(member,histories_team),tally);
        }
    }

  private:

    int d_league_size;
    const int d_max_history_length;
    const int d_num_histories;
    const random_scalar_view   d_start_cdf;
    const random_scalar_view   d_start_wt;
    const MC_Data_Texture_View d_mc_data;
    const scalar_view          d_y;
    const random_scalar_view   d_coeffs;
    random_scalar_view_2d      d_randoms;

};



} // end namespace alea

#endif // mc_solvers_MC_Components_hh

//---------------------------------------------------------------------------//
//                 end of MC_Components.hh
//---------------------------------------------------------------------------//
