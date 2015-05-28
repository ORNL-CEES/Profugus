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
    {}

    scalar_view weight;
    ord_view    state;
    ord_view    starting_ind;
    ord_view    row_length;
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
        }
    }

  private:

    const random_scalar_view   d_randoms;
    const History_Data         d_hist_data;
    const MC_Data_Texture_View d_mc_data;
};

//===========================================================================//
/*!
 * \class BinnedStateTransition
 * \brief Kernel for transitioning a new state
 *
 * This class is designed to be used within a Kokkos::parallel_for kernel
 * using a Kokkos::TeamPolicy.  This kernel differs from the previous one
 * in that it tries to achieve locality in data access by having each thread
 * team only process entries within a certain range of states.
 */
//===========================================================================//

class BinnedStateTransition
{
  public:

    typedef Kokkos::TeamPolicy<DEVICE> policy_type;
    typedef policy_type::member_type   policy_member;

    unsigned team_shmem_size(int team_size) const
    {
        return shared_ord_view::shmem_size(team_size);
    }

    BinnedStateTransition(scalar_view         randoms,
                          ord_view            old_states,
                          const History_Data &hist_data,
                          const MC_Data_View &mc_data,
                          int N)
        : d_randoms(randoms)
        , d_old_states(old_states)
        , d_states(hist_data.state)
        , d_weights(hist_data.weight)
        , d_P(mc_data.P)
        , d_W(mc_data.W)
        , d_inds(mc_data.inds)
        , d_offsets(mc_data.offsets)
        , d_num_histories(randoms.size())
        , d_N(N)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const policy_member &member) const
    {
        // Determine range of states to be processed by this team
        int state_begin, state_end;
        {
            int states_per_team = d_N / member.league_size();
            int extra_states = d_N % member.league_size();
            state_begin = member.league_rank() * states_per_team;
            if( member.league_rank() < extra_states )
            {
                state_begin += member.league_rank();
                states_per_team++;
            }
            else
            {
                state_begin += extra_states;
            }
            state_end = state_begin + states_per_team;
        }

        // Shared space for storing ids to be processed by team
        shared_ord_view ready_ids(  member.team_shmem(), member.team_size());

        int count = 0;
        int tid, hist, myindex, this_count;

        //
        // Loop over ALL histories, process those in specified range
        //

        // Process in batches of the team size
        for( int i=0; i<d_num_histories; i+=member.team_size() )
        {
            tid = i + member.team_rank();
            hist = -1;
            if( tid < d_num_histories        &&
                d_old_states(tid) >= state_begin &&
                d_old_states(tid) <  state_end )
            {
                hist = tid;
            }

            myindex = member.team_scan(static_cast<int>(hist>=0));
            this_count = myindex + (hist>=0);
            member.team_broadcast(this_count,member.team_size()-1);

            // Add my work to ready list if it is a valid index
            if( hist>=0 && count+myindex < member.team_size() )
                ready_ids(count+myindex) = tid;

            // If enough histories are ready, process them
            if( count+this_count >= member.team_size() )
            {
                // Process available histories
                member.team_barrier();
                process_history(ready_ids(member.team_rank()));

                // Populate ready list from work that didn't fit into list
                if( hist >= 0 )
                {
                    ready_ids((count+myindex)%member.team_size())=tid;
                }
                count = (count+this_count)%member.team_size();

            }
            else
            {
                count += this_count;
            }
        }

        // Process any remaining histories
        if( member.team_rank() < count )
        {
            process_history(ready_ids(member.team_rank()));
        }
    }

  private:

    KOKKOS_INLINE_FUNCTION
    void process_history(int tid) const
    {
        // Perform lower_bound search to get new state
        int start_ind  = d_offsets(d_states(tid));
        int row_length = d_offsets(d_states(tid)+1)-start_ind;
        int new_ind = lower_bound(d_P,start_ind,row_length,d_randoms(tid));

        // Update state
        if( new_ind == start_ind + row_length )
        {
            d_weights(tid) = 0.0;
            d_states(tid)  = -1;
        }
        else
        {
            d_weights(tid) *= d_W(new_ind);
            d_states(tid)   = d_inds(new_ind);
        }
    }


    const random_scalar_view   d_randoms;
    const random_ord_view      d_old_states;
    const ord_view             d_states;
    const scalar_view          d_weights;
    const random_scalar_view   d_P;
    const random_scalar_view   d_W;
    const random_ord_view      d_inds;
    const random_ord_view      d_offsets;
    int d_num_histories;
    int d_N;
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
    }

  private:

    const History_Data       d_hist_data;
    const random_scalar_view d_coeffs;
    const scalar_view        d_y;
};

//===========================================================================//
/*!
 * \class SharedMemTransition
 * \brief Kernel for transitioning a new state
 *
 * This class is designed to be used within a Kokkos::parallel_for kernel
 * using a Kokkos::TeamPolicy.
 */
//===========================================================================//

class SharedMemTransition
{
  public:

    typedef Kokkos::TeamPolicy<DEVICE> policy_type;
    typedef policy_type::member_type   policy_member;

    unsigned team_shmem_size(int team_size) const
    {
        return shared_ord_view::shmem_size(team_size) +
               shared_scalar_view::shmem_size(d_max_block_size);
    }

    int num_blocks() const{ return d_num_blocks; }

    SharedMemTransition(scalar_view         randoms,
                        ord_view            old_states,
                        const History_Data &hist_data,
                        const MC_Data_View &mc_data,
                        int N, int num_shared_values)
        : d_randoms(randoms)
        , d_old_states(old_states)
        , d_states(hist_data.state)
        , d_weights(hist_data.weight)
        , d_P(mc_data.P)
        , d_W(mc_data.W)
        , d_inds(mc_data.inds)
        , d_offsets(mc_data.offsets)
        , d_num_histories(randoms.size())
        , d_N(N)
    {
        // Copy offsets array to host
        REQUIRE( d_offsets.size() == (d_N+1) );
        ord_host_mirror offsets_host = Kokkos::create_mirror_view(d_offsets);
        Kokkos::deep_copy(offsets_host,d_offsets);

        // Determine how many "blocks" we need so that each block of P
        // has at most "num_shared_values" entries
        std::vector<int> states_per_block(1,0);
        states_per_block.reserve(d_N/num_shared_values);
        int sum = 0;
        d_max_block_size = 0;
        for( int i=0; i<d_N; ++i )
        {
            int num_this_row = offsets_host(i+1) - offsets_host(i);
            REQUIRE( num_this_row <= num_shared_values );
            if( sum + num_this_row <= num_shared_values )
            {
                states_per_block.back()++;
                sum += num_this_row;
            }
            else
            {
                d_max_block_size =
                    sum>d_max_block_size ? sum : d_max_block_size;
                states_per_block.push_back(1);
                sum = num_this_row;
            }
        }
        ENSURE( d_max_block_size <= num_shared_values );
        std::cout << "Maximum block size is " << d_max_block_size <<
            ", max size allowed is " << num_shared_values << std::endl;

        d_num_blocks = states_per_block.size();
        Kokkos::resize(d_block_offsets,d_num_blocks+1);
        ord_host_mirror block_offsets_host =
            Kokkos::create_mirror_view(d_block_offsets);
        block_offsets_host(0) = 0;
        for( int i=0; i<d_num_blocks; ++i )
        {
            REQUIRE( states_per_block[i] <= num_shared_values );
            block_offsets_host(i+1) = block_offsets_host(i) +
                states_per_block[i];
        }

        Kokkos::deep_copy(d_block_offsets,block_offsets_host);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const policy_member &member) const
    {
        shared_scalar_view P_shared(member.team_shmem(),d_max_block_size);

        // Determine range of states to be processed by this team
        int state_begin = d_block_offsets(member.league_rank());
        int state_end   = d_block_offsets(member.league_rank()+1);
        int block_offset = d_offsets(state_begin);
        int block_end = d_offsets(state_end);
        member.team_barrier();
        {
            for( int i=block_offset; i<block_end; i+=member.team_size() )
            {
                if( i + member.team_rank() < block_end )
                {
                    P_shared(member.team_rank()+i-block_offset) =
                        d_P(i+member.team_rank());
                }
            }
        }

        // Shared space for storing ids to be processed by team
        shared_ord_view ready_ids(member.team_shmem(), member.team_size());

        int count = 0;
        int tid, hist, myindex, this_count;

        //
        // Loop over ALL histories, process those in specified range
        //

        // Process in batches of the team size
        for( int i=0; i<d_num_histories; i+=member.team_size() )
        {
            tid = i + member.team_rank();
            hist = -1;
            if( tid < d_num_histories        &&
                d_old_states(tid) >= state_begin &&
                d_old_states(tid) <  state_end )
            {
                hist = tid;
            }

            myindex = member.team_scan(static_cast<int>(hist>=0));
            this_count = myindex + (hist>=0);
            member.team_broadcast(this_count,member.team_size()-1);

            // Add my work to ready list if it is a valid index
            if( hist>=0 && count+myindex < member.team_size() )
                ready_ids(count+myindex) = tid;

            // If enough histories are ready, process them
            if( count+this_count >= member.team_size() )
            {
                // Process available histories
                member.team_barrier();
                process_history(ready_ids(member.team_rank()),P_shared,
                    block_offset);

                // Populate ready list from work that didn't fit into list
                if( hist >= 0 )
                {
                    ready_ids((count+myindex)%member.team_size())=tid;
                }
                count = (count+this_count)%member.team_size();

            }
            else
            {
                count += this_count;
            }
        }

        // Process any remaining histories
        if( member.team_rank() < count )
        {
            member.team_barrier();
            process_history(ready_ids(member.team_rank()),P_shared,
                block_offset);
        }
    }

  private:

    KOKKOS_INLINE_FUNCTION
    void process_history(int tid, shared_scalar_view &P_shared,
            int block_offset) const
    {
        // Starting index and length of current row
        int start_ind  = d_offsets(d_states(tid));
        int row_length = d_offsets(d_states(tid)+1)-start_ind;

        // Perform lower_bound search to get new state
        // Decrease starting index by block_offset to make the index relative
        // to local block (P_shared is only a small portion of the global P)
        int new_ind = lower_bound(P_shared,start_ind-block_offset,row_length,
            d_randoms(tid));

        // Put index back into global index space
        new_ind += block_offset;

        // Update state
        if( new_ind == start_ind + row_length )
        {
            d_weights(tid) = 0.0;
            d_states(tid)  = -1;
        }
        else
        {
            d_weights(tid) *= d_W(new_ind);
            d_states(tid)   = d_inds(new_ind);
        }
    }

    const random_scalar_view   d_randoms;
    const random_ord_view      d_old_states;
    const ord_view             d_states;
    const scalar_view          d_weights;
    ord_view                   d_block_offsets;
    const random_scalar_view   d_P;
    const random_scalar_view   d_W;
    const random_ord_view      d_inds;
    const random_ord_view      d_offsets;
    int d_num_histories;
    int d_N;
    int d_num_blocks;
    int d_max_block_size;
};



} // end namespace alea

#endif // mc_solvers_MC_Components_hh

//---------------------------------------------------------------------------//
//                 end of MC_Components.hh
//---------------------------------------------------------------------------//
