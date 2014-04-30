//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Rebalance.cc
 * \author Thomas M. Evans
 * \date   Thu Feb 28 13:05:02 2013
 * \brief  Continuous_Fission_Rebalance template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Rebalance_cc
#define mc_Fission_Rebalance_cc

#include "Continuous_Fission_Rebalance.hh"

#include <sstream>
#include <numeric>
#include "harness/DBC.hh"
#include "comm/Timing.hh"

namespace shift
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Physics>
Continuous_Fission_Rebalance<Physics>::Continuous_Fission_Rebalance(
    SP_Boundary_Mesh bnd_mesh)
    : d_bnd_mesh(bnd_mesh)
    , d_num_sets(d_bnd_mesh->num_sets())
    , d_set(d_bnd_mesh->set())
    , d_num_blocks(d_bnd_mesh->num_blocks())
    , d_block(d_bnd_mesh->block())
    , d_left(-1)
    , d_right(-1)
    , d_num_nbors(0)
    , d_target(0)
    , d_target_set(0)
    , d_size_fs(Physics::fission_site_bytes())
{
    Require (d_bnd_mesh);

    // return if we are on 1 set
    if (d_num_sets == 1)
        return;

    // allocate space for the number of fission sites across all sets
    d_sites_set.resize(d_num_sets);

    // allocate space for the number of fission sites for all blocks on each
    // set
    d_sites_block.resize(d_num_blocks);
    d_block_surplus.resize(d_num_blocks);

    // >>> BLOCK-CONSTANT intra-set communicator
    //     make communicator for all sets that have the same block, ie. block
    //     0 in sets 1,...,N and block 1 in sets 1,...,N each have their own
    //     communicator
    nemesis::split(d_block, 0, d_block_const_comm);

    // >>> SET-CONSTANT communicator
    //     make communicator for all blocks within a set, ie. blocks 0,N in
    //     set 1 has a communicator, blocks 0,N in set 2 has a communicator,
    //     etc.
    nemesis::split(d_set, 0, d_set_const_comm);

    // setup left/right neighbor communicators
    nemesis::set_internal_comm(d_block_const_comm);
    Check (nemesis::node()  == d_set);
    Check (nemesis::nodes() == d_num_sets);

    // get left and right neighbors
    if (d_set > 0)
    {
        d_left = d_set - 1;
        ++d_num_nbors;
    }
    if (d_set < d_num_sets - 1)
    {
        d_right = d_set + 1;
        ++d_num_nbors;
    }
    Check (d_num_nbors > 0);

    nemesis::reset_internal_comm();
    Ensure (nemesis::nodes() == d_bnd_mesh->num_domains());
    Ensure (nemesis::node() == d_bnd_mesh->domain());
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Rebalance the fission bank across sets.
 *
 * When the number of sets is greater than 1, the fission bank is rebalanced
 * across all of the sets using the algorithm described in the
 * Continuous_Fission_Rebalance class description.
 *
 * When the number of sets is equal to 1,
 *
 * \param fission_bank on entry the fission bank on this set/block; on exit a
 * rebalanced fission bank
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::rebalance(
    Fission_Site_Container_t &fission_bank)
{
    Require (d_bnd_mesh);
    Require (d_bnd_mesh->num_sets() == d_num_sets);
    Require (d_bnd_mesh->set() == d_set);

    SCOPED_TIMER("shift::Continuous_Fission_Rebalance.rebalance");

    // return if only 1 set
    if (d_num_sets == 1)
    {
        // the target number for this set is simply the number of sites summed
        // over all the blocks
        d_target_set = fission_bank.size();
        d_target     = fission_bank.size();

        nemesis::global_sum(d_target_set);

        // number of global fissions is the same as the number on the set
        d_num_global = d_target_set;

        // we are done
        return;
    }

    // >>> CHANGING COMMUNICATORS TO BLOCK-CONST COMMUNICATION
    nemesis::set_internal_comm(d_block_const_comm);
    Check (nemesis::node()  == d_set);
    Check (nemesis::nodes() == d_num_sets);

    // set-up global/local fission bank parameters
    fission_bank_parameters(fission_bank);
    nemesis::global_barrier();

    // initialize send/receive counters for this rebalance
    d_num_recv = 0;
    d_num_send = 0;
    d_num_iter = 0;

    // actual fissions on the set
    int set_sites = 0;

    // iterate until the fission bank is balanced
    int converged_sets = 0;
    while (converged_sets < d_num_sets)
    {
        // send/receive particles for this iteration
        communicate(fission_bank);

        // calculate the number of fissions on the set
        set_sites = fission_bank.size();
        nemesis::set_internal_comm(d_set_const_comm);
        nemesis::global_sum(&set_sites, 1);

        // check for convergence
        nemesis::set_internal_comm(d_block_const_comm);
        converged_sets = 0;
        if (set_sites == d_target_set) converged_sets = 1;
        nemesis::global_sum(converged_sets);

        // reset internal comm
        nemesis::reset_internal_comm();

        // update the number of sites on each domain if we are not converged
        if (converged_sets < d_num_sets * d_num_blocks)
        {
            calc_num_sites(fission_bank);
        }

        // iteration counter
        ++d_num_iter;
    }

    nemesis::reset_internal_comm();
    // <<< COMMUNICATION RESET TO NORMAL COMMUNICATOR

    // set the target on each domain (number of fission sites after rebalance)
    d_target = fission_bank.size();

    Ensure (nemesis::nodes() == d_bnd_mesh->num_domains());
    Ensure (nemesis::node() == d_bnd_mesh->domain());
    Ensure (d_num_blocks == 1 ? d_target_set == d_target :
            d_target <= d_target_set);

#ifdef ENSURE_ON
    int global_check = fission_bank.size();
    nemesis::global_sum(global_check);
    Validate(global_check == d_num_global,
          "Failed to preserve global fission sites: Calculated = "
          << global_check << "; Expected = " << d_num_global
          << "; Set = " << d_set << "; Block = " << d_block
          << "Set target = " << d_target_set
          << "; Actual on block = " << fission_bank.size());
#endif
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the global/local parameters of the fission bank.
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::fission_bank_parameters(
    const Fission_Site_Container_t &fission_bank)
{
    Require (d_num_sets > 1);
    Require (d_sites_set.size() == d_num_sets);

    // calculate the number of fission sites and local array bounds
    calc_num_sites(fission_bank);

    // determine the global number of fission sites
    d_num_global = std::accumulate(d_sites_set.begin(), d_sites_set.end(), 0);
    Check (d_num_global > 0);

    // calculate the target number of fission sites on this set
    d_target_set = d_num_global / d_num_sets;

    // initialize the target low-edge (first) array boundary on this set
    d_target_bnds.first = d_target_set * d_set;

    // determine extra sites when for non-uniform numbers of sites
    int pad = d_num_global - d_target_set * d_num_sets;
    Check (pad >= 0 && pad < d_num_sets);

    // add sites to account for padding, one site is added to each set until
    // the correct global number of sites is attained
    if (d_set < pad)
    {
        ++d_target_set;
        d_target_bnds.first += d_set;
    }
    else
    {
        d_target_bnds.first += pad;
    }

    // calculate the high-edge (last) array boundary on this set
    d_target_bnds.second = d_target_bnds.first + d_target_set - 1;

#ifdef CHECK_ON
    nemesis::set_internal_comm(d_block_const_comm);
    int global_check = d_target_set;
    nemesis::global_sum(global_check);
    nemesis::reset_internal_comm();
    Validate(global_check == d_num_global,
            "Failed to accurately pad sets for non-uniform fission sites.");
#endif

    Ensure (d_target_bnds.second - d_target_bnds.first + 1 == d_target_set);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Communicate fission bank sites during a rebalance step.
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::communicate(
    Fission_Site_Container_t &fission_bank)
{
    Require (d_recv_left.empty());
    Require (d_recv_right.empty());

    // send/receive from right/left
    int num_send_left  = 0;
    int num_send_right = 0;
    int num_recv_left  = 0;
    int num_recv_right = 0;

    // calculate the number of sites to send/recv; sends are constrained by
    // the number of sites currently on the set

    // determine the number to send to the left
    if (d_bnds.first < d_target_bnds.first)
    {
        Check (d_left != -1);
        num_send_left = std::min(d_target_bnds.first - d_bnds.first,
                                 d_sites_set[d_set]);
    }
    else if (d_bnds.first > d_target_bnds.first)
    {
        Check (d_left != -1);
        num_recv_left = std::min(d_bnds.first - d_target_bnds.first,
                                 d_sites_set[d_set - 1]);
    }

    // determine the number to send to the right
    if (d_bnds.second > d_target_bnds.second)
    {
        Check (d_right != -1);
        num_send_right = std::min(d_bnds.second - d_target_bnds.second,
                                  d_sites_set[d_set]);
    }

    // determine the number we receive from the right
    else if (d_bnds.second < d_target_bnds.second)
    {
        Check (d_right != -1);
        num_recv_right = std::min(d_target_bnds.second - d_bnds.second,
                                  d_sites_set[d_set + 1]);
    }

    Check (num_send_left >= 0 && num_send_left <= d_sites_set[d_set]);
    Check (num_send_right >= 0 && num_send_right <= d_sites_set[d_set]);
    Check (num_send_left + num_send_right <= d_sites_set[d_set]);
    Check (num_send_left ? num_recv_left == 0 : true);
    Check (num_send_right ? num_recv_right == 0 : true);

    // set the number of sites sent left and right on the full set
    int num_send_left_set  = num_send_left;
    int num_send_right_set = num_send_right;

    // rebalance across the blocks --> num_send after this call now refers to
    // the number of sites sent on each block within the set
    block_rebalance(num_send_left);
    block_rebalance(num_send_right);

    // set to block-to-block communication across sets
    nemesis::set_internal_comm(d_block_const_comm);

    // determine the number to receive on each block
    determine_block_rcvs(num_send_left, num_send_right,
                         num_send_left_set, num_send_right_set,
                         num_recv_left, num_recv_right);

    // post receives
    post_receives(num_recv_left, d_recv_left, d_left, d_handle_left, 303);
    post_receives(num_recv_right, d_recv_right, d_right, d_handle_right, 304);

    // send sites
    send(num_send_left, fission_bank, d_left, 304);
    send(num_send_right, fission_bank, d_right, 303);

    // receive sites
    receive(num_recv_left, fission_bank, d_recv_left, d_left, d_handle_left,
            303);
    receive(num_recv_right, fission_bank, d_recv_right, d_right,
            d_handle_right, 304);

    Ensure (!d_handle_right.inuse());
    Ensure (!d_handle_left.inuse());

    // reset to internal comm
    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of physics sites on each set and the current
 * fission bank array bounds.
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::calc_num_sites(
    const Fission_Site_Container_t &fission_bank)
{
    Require (d_sites_set.size() == d_num_sets);
    Require (d_sites_block.size() == d_num_blocks);

    // do a global reduction to determine the number of sites on this set
    nemesis::set_internal_comm(d_set_const_comm);
    Check (nemesis::nodes() == d_num_blocks);
    std::fill(d_sites_block.begin(), d_sites_block.end(), 0);
    d_sites_block[d_block] = fission_bank.size();
    nemesis::global_sum(&d_sites_block[0], d_num_blocks);

    // do a global reduction to determine the current number of sites on all
    // sets
    nemesis::set_internal_comm(d_block_const_comm);
    Check (nemesis::nodes() == d_num_sets);
    std::fill(d_sites_set.begin(), d_sites_set.end(), 0);
    d_sites_set[d_set] = std::accumulate(d_sites_block.begin(),
                                         d_sites_block.end(), 0);
    nemesis::global_sum(&d_sites_set[0], d_num_sets);

    // reset the communicator
    nemesis::reset_internal_comm();

    // make the array bounds on this set --> the array bounds are (first,last)
    d_bnds.first  = std::accumulate(&d_sites_set[0], &d_sites_set[0]+d_set, 0);
    d_bnds.second = d_bnds.first + d_sites_set[d_set] - 1;
    Check (d_bnds.second - d_bnds.first + 1 == d_sites_set[d_set]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Rebalance across blocks in a set.
 *
 * \param num_send on input it is the number to send over the whole set; on
 * output it is the number to send from this block
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::block_rebalance(int &num_send)
{
    Require (d_sites_block.size() == d_block_surplus.size());

    // return immediately if there is only a single block per set
    if (d_num_blocks == 1)
        return;

    // return if there are no sites to send
    if (num_send == 0)
        return;

    // all of the following calculations should be per set
    nemesis::set_internal_comm(d_set_const_comm);

    // calculate the surplus sites per block (deviation from average)
    std::copy(d_sites_block.begin(), d_sites_block.end(),
              d_block_surplus.begin());
    double avg = static_cast<double>(
        std::accumulate(d_sites_block.begin(), d_sites_block.end(), 0)) /
                 d_num_blocks;
    for (int b = 0; b < d_num_blocks; ++b)
    {
        d_block_surplus[b] = std::max(0.0, d_block_surplus[b] - avg);
    }

    // the number of sites above the average number per block
    int num_surplus = std::accumulate(d_block_surplus.begin(),
                                      d_block_surplus.end(), 0);

    // running counter for sends
    int cntr = 0;

    // if the number of surplus is large enough to accumodate the send size,
    // remove sites from each surplus block one at a time until we have enough
    if (num_surplus > num_send)
    {
        // initialize the running counter
        cntr = num_surplus;

        // iterate through the blocks until we get the correct number of sites
        // to send
        while (cntr > num_send)
        {
            for (int b = 0; b < d_num_blocks; ++b)
            {
                // only pull from positive surplus
                if (d_block_surplus[b] > 0)
                {
                    --d_block_surplus[b];
                    --cntr;
                }

                // break out of loop when we have enough sites
                if (cntr == num_send)
                    break;
            }
        }
    }

    // if the number of surplus is not large enough we have to pull from other
    // blocks
    if (num_surplus < num_send)
    {
        // initialize the running counter
        cntr = num_surplus;

        // iterate thorugh teh blocks until we get the correct number of sites
        // to send
        while (cntr < num_send)
        {
            for (int b = 0; b < d_num_blocks; ++b)
            {
                // only pull from blocks that have still have sites
                if (d_sites_block[b] - d_block_surplus[b] > 0)
                {
                    ++d_block_surplus[b];
                    ++cntr;
                }

                // break out of loop when we have enough sites
                if (cntr == num_send)
                    break;
            }
        }
    }

    Check (std::accumulate(d_block_surplus.begin(), d_block_surplus.end(), 0)
           == num_send);

    // update the number of sites per block by subracting the number that will
    // be sent
    for (int b = 0; b < d_num_blocks; ++b)
    {
        d_sites_block[b] -= d_block_surplus[b];
        Check (d_sites_block[b] >= 0);
    }

    // reset the number to send on this block
    num_send = d_block_surplus[d_block];

    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the number of sites to receive on this block.
 *
 * \param num_recv_left on input it is the number of sites to receive on the
 * set from the neighbor to the left; on output it is the number of sites that
 * will be received on this block
 *
 * \param num_recv_left on input it is the number of sites to receive on the
 * set from the neighbor to the right; on output it is the number of sites
 * that will be received on this block
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::determine_block_rcvs(
    int  num_send_left,
    int  num_send_right,
    int  num_send_left_set,
    int  num_send_right_set,
    int &num_recv_left,
    int &num_recv_right)
{
    Require (nemesis::nodes() == d_num_sets);
    Require (num_send_left <= num_send_left_set);
    Require (num_send_right <= num_send_right_set);
    Require (!d_handle_left.inuse());
    Require (!d_handle_right.inuse());

    // return if there is only 1 block
    if (d_num_blocks == 1)
    {
        Check (num_send_left == num_send_left_set);
        Check (num_send_right == num_send_right_set);
        return;
    }

    // local store of number received on the set
    int num_recv_left_set  = num_recv_left;
    int num_recv_right_set = num_recv_right;

    // if any block on this set receives from the right/left then post a
    // receive on this block (the size received could be zero)
    if (num_recv_left_set)
    {
        // post a receive from the left neighbor
        nemesis::receive_async(
            d_handle_left, &num_recv_left, 1, d_left, 401);
    }
    if (num_recv_right_set)
    {
        // post a receive from the left neighbor
        nemesis::receive_async(
            d_handle_right, &num_recv_right, 1, d_right, 402);
    }

    // send to the right/left
    if (num_send_left_set)
    {
        nemesis::send(&num_send_left, 1, d_left, 402);
    }
    if (num_send_right_set)
    {
        nemesis::send(&num_send_right, 1, d_right, 401);
    }

    // wait on receives (if the handle is unassigned this is a no-op)
    d_handle_left.wait();
    d_handle_right.wait();


    Ensure (num_recv_left >= 0);
    Ensure (num_recv_right >= 0);
    Ensure (!d_handle_left.inuse());
    Ensure (!d_handle_right.inuse());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post receives.
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::post_receives(
    int                       num_recv,
    Fission_Site_Container_t &recv_bank,
    int                       destination,
    nemesis::Request         &handle,
    int                       tag)
{
    if (num_recv)
    {
        Require (destination >= 0 && destination < d_num_sets);

        // allocate space in the receive buffer
        recv_bank.resize(num_recv);

        // make a void * to the container
        void *buffer = &recv_bank[0];

        // post the receive
        nemesis::receive_async(handle, reinterpret_cast<char *>(buffer),
                               num_recv * d_size_fs, destination, tag);

        // increment receive counter
        ++d_num_recv;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send sites.
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::send(
    int                       num_send,
    Fission_Site_Container_t &bank,
    int                       destination,
    int                       tag)
{
    Require (bank.size() >= num_send);
    Remember (int size = bank.size());

    if (num_send)
    {
        Require (!bank.empty());

        // make a void * to the LAST num_send sites in the bank
        const void *buffer = &bank[0] + (bank.size() - num_send);

        // send the first *num_send* sites in the fission bank
        nemesis::Request handle = nemesis::send_async(
            reinterpret_cast<const char *>(buffer), num_send * d_size_fs,
            destination, tag);
        handle.wait();

        // pop sites off the bank that have been sent
        for (int n = 0; n < num_send; ++n)
        {
            Check (!bank.empty());
            bank.pop_back();
        }

        // increment send counter
        ++d_num_send;
    }

    Ensure (bank.size() == size - num_send);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Receive sites.
 */
template<class Physics>
void Continuous_Fission_Rebalance<Physics>::receive(
    int                       num_recv,
    Fission_Site_Container_t &bank,
    Fission_Site_Container_t &recv_bank,
    int                       destination,
    nemesis::Request         &handle,
    int                       tag)
{
    Remember (int size = bank.size());

   if (num_recv)
    {
        Require (destination >= 0 && destination < d_num_sets);
        Require (handle.inuse());

        // wait for the data to arrive
        handle.wait();

        // append sites to the bank
        for (fs_iterator s = recv_bank.begin(); s != recv_bank.end(); ++s)
        {
            bank.push_back(*s);
        }

        // clear the receive buffer
        recv_bank.clear();
    }

    Ensure (!handle.inuse());
    Ensure (bank.size() == size + num_recv);
}

} // end namespace shift

#endif // mc_Fission_Rebalance_cc

//---------------------------------------------------------------------------//
//                 end of Fission_Rebalance.cc
//---------------------------------------------------------------------------//
