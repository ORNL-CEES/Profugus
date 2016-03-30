//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Rebalance.t.hh
 * \author Thomas M. Evans
 * \date   Monday May 5 11:37:31 2014
 * \brief  Fission_Rebalance template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fission_Rebalance_t_cuh
#define cuda_mc_Fission_Rebalance_t_cuh

#include "Fission_Rebalance.hh"

#include <algorithm>
#include <sstream>
#include <numeric>

#include "cuda_utils/CudaDBC.hh"
#include "comm/Timing.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Fission_Rebalance<Geometry>::Fission_Rebalance()
    : d_num_sets(profugus::nodes())
    , d_set(profugus::node())
    , d_left(-1)
    , d_right(-1)
    , d_num_nbors(0)
    , d_target_set(0)
    , d_sites_set(d_num_sets)
    , d_size_fs(Physics_t::fission_site_bytes())
{
    // return if we are on 1 set
    if (d_num_sets == 1)
        return;

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
    CHECK(d_num_nbors > 0);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Rebalance the fission bank across sets.
 *
 * When the number of sets is greater than 1, the fission bank is rebalanced
 * across all of the sets using the algorithm described in the
 * Fission_Rebalance class description.
 *
 * When the number of sets is equal to 1, this is a no-op.
 *
 * \param fission_bank on entry the fission bank on this set/block; on exit a
 * rebalanced fission bank
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::rebalance(
        Fission_Site_Container_t &fission_bank)
{
    SCOPED_TIMER("CUDA_MC::Fission_Rebalance.rebalance");

    // initialize send/receive counters for this rebalance
    d_num_recv = 0;
    d_num_send = 0;
    d_num_iter = 0;

    // return if only 1 set
    if (d_num_sets == 1)
    {
        // the target number for this set is simply the number of sites
        d_target_set = fission_bank.size();

        // number of global fissions is the same as the number on the set
        d_num_global = d_target_set;

        // we are done
        return;
    }

    // set-up global/local fission bank parameters
    fission_bank_parameters(fission_bank);
    profugus::global_barrier();

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

        // check for convergence
        converged_sets = 0;
        if (set_sites == d_target_set) converged_sets = 1;
        profugus::global_sum(converged_sets);

        // update the number of sites on each domain if we are not converged
        if (converged_sets < d_num_sets)
        {
            calc_num_sites(fission_bank);
        }

        // iteration counter
        ++d_num_iter;
    }

    // set the target on each domain (number of fission sites after rebalance)
    d_target_set = fission_bank.size();

#ifdef ENSURE_ON
    int global_check = fission_bank.size();
    profugus::global_sum(global_check);
    VALIDATE(global_check == d_num_global,
             "Failed to preserve global fission sites: Calculated = "
             << global_check << "; Expected = " << d_num_global
             << "; Set = " << d_set
             << "; Set target = " << d_target_set
             << "; Actual on set = " << fission_bank.size());
#endif
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the global/local parameters of the fission bank.
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::fission_bank_parameters(
    const Fission_Site_Container_t &fission_bank)
{
    REQUIRE(d_num_sets > 1);
    REQUIRE(d_sites_set.size() == d_num_sets);

    // calculate the number of fission sites and local array bounds
    calc_num_sites(fission_bank);

    // determine the global number of fission sites
    d_num_global = std::accumulate(d_sites_set.begin(), d_sites_set.end(), 0);
    CHECK(d_num_global > 0);

    // calculate the target number of fission sites on this set
    d_target_set = d_num_global / d_num_sets;

    // initialize the target low-edge (first) array boundary on this set
    d_target_bnds.first = d_target_set * d_set;

    // determine extra sites when for non-uniform numbers of sites
    int pad = d_num_global - d_target_set * d_num_sets;
    CHECK(pad >= 0 && pad < d_num_sets);

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
    int global_check = d_target_set;
    profugus::global_sum(global_check);
    VALIDATE(global_check == d_num_global,
            "Failed to accurately pad sets for non-uniform fission sites.");
#endif

    ENSURE(d_target_bnds.second - d_target_bnds.first + 1 == d_target_set);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Communicate fission bank sites during a rebalance step.
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::communicate(
        Fission_Site_Container_t &fission_bank)
{
    REQUIRE(d_recv_left.empty());
    REQUIRE(d_recv_right.empty());

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
        CHECK(d_left != -1);
        num_send_left = std::min(d_target_bnds.first - d_bnds.first,
                                 d_sites_set[d_set]);
    }
    else if (d_bnds.first > d_target_bnds.first)
    {
        CHECK(d_left != -1);
        num_recv_left = std::min(d_bnds.first - d_target_bnds.first,
                                 d_sites_set[d_set - 1]);
    }

    // determine the number to send to the right
    if (d_bnds.second > d_target_bnds.second)
    {
        CHECK(d_right != -1);
        num_send_right = std::min(d_bnds.second - d_target_bnds.second,
                                  d_sites_set[d_set]);
    }

    // determine the number we receive from the right
    else if (d_bnds.second < d_target_bnds.second)
    {
        CHECK(d_right != -1);
        num_recv_right = std::min(d_target_bnds.second - d_bnds.second,
                                  d_sites_set[d_set + 1]);
    }

    CHECK(num_send_left >= 0 && num_send_left <= d_sites_set[d_set]);
    CHECK(num_send_right >= 0 && num_send_right <= d_sites_set[d_set]);
    CHECK(num_send_left + num_send_right <= d_sites_set[d_set]);
    CHECK(num_send_left ? num_recv_left == 0 : true);
    CHECK(num_send_right ? num_recv_right == 0 : true);

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

    ENSURE(!d_handle_right.inuse());
    ENSURE(!d_handle_left.inuse());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of physics sites on each set and the current
 * fission bank array bounds.
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::calc_num_sites(
    const Fission_Site_Container_t &fission_bank)
{
    REQUIRE(d_sites_set.size() == d_num_sets);
    REQUIRE(profugus::nodes() == d_num_sets);

    // do a global reduction to determine the current number of sites on all
    // sets
    std::fill(d_sites_set.begin(), d_sites_set.end(), 0);
    d_sites_set[d_set] = fission_bank.size();
    profugus::global_sum(&d_sites_set[0], d_num_sets);

    // make the array bounds on this set --> the array bounds are (first,last)
    d_bnds.first  = std::accumulate(&d_sites_set[0], &d_sites_set[0]+d_set, 0);
    d_bnds.second = d_bnds.first + d_sites_set[d_set] - 1;
    CHECK(d_bnds.second - d_bnds.first + 1 == d_sites_set[d_set]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post receives.
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::post_receives(
        int                       num_recv,
        Fission_Site_Container_t &recv_bank,
        int                       destination,
        profugus::Request         &handle,
        int                       tag)
{
    if (num_recv)
    {
        REQUIRE(destination >= 0 && destination < d_num_sets);

        // allocate space in the receive buffer
        recv_bank.resize(num_recv);

        // make a void * to the container
        void *buffer = &recv_bank[0];

        // post the receive
        profugus::receive_async(handle, reinterpret_cast<char *>(buffer),
                               num_recv * d_size_fs, destination, tag);

        // increment receive counter
        ++d_num_recv;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send sites.
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::send(
        int                       num_send,
        Fission_Site_Container_t &bank,
        int                       destination,
        int                       tag)
{
    REQUIRE(bank.size() >= num_send);

    if (num_send)
    {
        REQUIRE(!bank.empty());

        // make a void * to the LAST num_send sites in the bank
        const void *buffer = &bank[0] + (bank.size() - num_send);

        // send the first *num_send* sites in the fission bank
        profugus::Request handle = profugus::send_async(
            reinterpret_cast<const char *>(buffer), num_send * d_size_fs,
            destination, tag);
        handle.wait();

        // pop sites off the bank that have been sent
        for (int n = 0; n < num_send; ++n)
        {
            CHECK(!bank.empty());
            bank.pop_back();
        }

        // increment send counter
        ++d_num_send;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Receive sites.
 */
template <class Geometry>
void Fission_Rebalance<Geometry>::receive(
        int                       num_recv,
        Fission_Site_Container_t &bank,
        Fission_Site_Container_t &recv_bank,
        int                       destination,
        profugus::Request        &handle,
        int                       tag)
{
   if (num_recv)
    {
        REQUIRE(destination >= 0 && destination < d_num_sets);
        REQUIRE(handle.inuse());

        // wait for the data to arrive
        handle.wait();

        // append sites to the bank
        for (auto s = recv_bank.begin(); s != recv_bank.end(); ++s)
        {
            bank.push_back(*s);
        }

        // clear the receive buffer
        recv_bank.clear();
    }

    ENSURE(!handle.inuse());
}

} // end namespace cuda_profugus

#endif // cuda_mc_Fission_Rebalance_t.hh

//---------------------------------------------------------------------------//
//                 end of Fission_Rebalance.t.hh
//---------------------------------------------------------------------------//
