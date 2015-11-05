//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Rebalance.hh
 * \author Thomas M. Evans
 * \date   Monday May 5 11:13:55 2014
 * \brief  Fission_Rebalance class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Rebalance_hh
#define mc_Fission_Rebalance_hh

#include <utility>

#include "comm/global.hh"
#include "utils/Definitions.hh"
#include "Physics.hh"

// Remove this once templated
#include "geometry/RTK_Geometry.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Rebalance
 * \brief Rebalance the fission bank across sets for continuous fission site
 * containers.
 *
 * The fission bank can go out of balance across multiple sets (each set has
 * the complete mesh, possibly decomposed into multiple blocks within each
 * set).  In order to keep the problem from becomming load-imbalanced, we need
 * to rebalance the fission bank.  Because each set has the complete problem
 * geometry, it does not matter which fission sites get moved between sets.
 *
 * This version of the fission rebalance algorithm \b assumes only 1-block per
 * set (ie. full domain replication).  The generalized Fission Bank algorithm
 * in Shift can handle multiple blocks per set. The algorithm we have
 * developed is a modification of Romano and Forgets fission bank algorithm
 * (Romano, P K, and B Forget. “Parallel Fission Bank Algorithms in Monte
 * Carlo Criticality Calculations.” \e Nuclear \e Science \e and \e
 * Engineering \b 170 125–135.).
 *
 * Consider a fission bank with 1003 sampled sites on 4 processors. The
 * distribution on each processor is given by \c n_i.  If we consider the
 * fission bank as a \e global array, on each processor the first/last bounds
 * are given by \c a_i/b_i where \c i is the processor id.  The \e balanced
 * array bounds are given by \c ax_i/bx_i.  The current state of the fission
 * bank is as follows:
 *
 * \verbatim
            i        a_i        b_i       ax_i       bx_i        n_i
   -----------------------------------------------------------------
            0          0        558          0        250        559
            1        559        672        251        501        114
            2        673        822        502        752        150
            3        823       1002        753       1002        180
   Total N =  1003
   \endverbatim
 *
 * Here, processors 0, 1, 2, and 3 have 559, 114, 150, and 180 sites.  The
 * desired distribution should be 251, 251, 251, and 250, respectively. We
 * need to rebalance, so in the first iteration we do the following:
 *
 * \verbatim
   Iteration =  0
   On 0, sending   308 particles to   1
   On 1, receiving 308 particles from 0
   On 1, sending   114 particles to   2
   On 2, receiving 114 particles from 1
   On 2, sending   70 particles to   3
   On 3, receiving 70 particles from 2

            i        a_i        b_i       ax_i       bx_i        n_i
   -----------------------------------------------------------------
            0          0        250          0        250        251
            1        251        558        251        501        308
            2        559        752        502        752        194
            3        753       1002        753       1002        250
   Total N =  1003
   \endverbatim
 *
 * At most, each set only performs 2 communications with its adjacent sets.
 * However, because this example was significantly out of balance, more
 * iterations are needed:
 *
 * \verbatim
   Iteration =  1
   On 1, sending   57 particles to   2
   On 2, receiving 57 particles from 1

            i        a_i        b_i       ax_i       bx_i        n_i
   -----------------------------------------------------------------
            0          0        250          0        250        251
            1        251        501        251        501        251
            2        502        752        502        752        251
            3        753       1002        753       1002        250
   Total N =  1003

   Finished in 2 iterations
   \endverbatim
 *
 * After, 2 iterations the problem is rebalanced.  The total number of
 * send/receives by set are:
 *
 * \verbatim
            i      Sends   Receives
   --------------------------------
            0          1          0
            1          2          1
            2          1          2
            3          0          1
   \endverbatim
 *
 * To summarize, in each iteration a set has at most 2 communications with its
 * nearest set neighbor.
 */
/*!
 * \example mc/test/tstFission_Rebalance.cc
 *
 * Test of Fission_Rebalance.
 */
//===========================================================================//

class Fission_Rebalance
{
  public:
    //@{
    //! Typedefs.
    typedef Physics<Core>                               Physics_t;
    typedef typename Physics_t::Fission_Site            Fission_Site_t;
    typedef typename Physics_t::Fission_Site_Container  Fission_Site_Container_t;
    typedef typename Physics_t::Space_Vector            Space_Vector;
    typedef def::Vec_Int                                Vec_Int;
    typedef std::pair<int, int>                         Array_Bnds;
    //@}

  public:
    // Constructor.
    Fission_Rebalance();

    // Rebalance the fission bank across all sets.
    void rebalance(Fission_Site_Container_t &fission_bank);

    // >>> ACCESSORS

    //! Number of fissions on this domain after a rebalance.
    int num_fissions() const { return d_target_set; }

    //! Total, global number of fissions (same before/after rebalance).
    int num_global_fissions() const { return d_num_global; }

    //! Get target bounds on the current set for this rebalance step.
    const Array_Bnds& target_array_bnds() const { return d_target_bnds; }

    //! Return number of sends during a rebalance step on this set.
    int num_sends() const { return d_num_send; }

    //! Return number of receives during a rebalance step on this set.
    int num_receives() const { return d_num_recv; }

    //! Return number of iterations for this rebalance.
    int num_iterations() const { return d_num_iter; }

  private:
    // >>> IMPLEMENTATION

    // Calculate global/local fission bank parameters.
    void fission_bank_parameters(const Fission_Site_Container_t &fission_bank);

    // Communicate fission bank sites during a rebalance step.
    void communicate(Fission_Site_Container_t &fission_bank);

    // Calculate the number of fission sites across all sets.
    void calc_num_sites(const Fission_Site_Container_t &fission_bank);

    // Post receives.
    void post_receives(int num_recv, Fission_Site_Container_t &recv_bank,
                       int destination, profugus::Request &handle, int tag);

    // Send.
    void send(int num_send, Fission_Site_Container_t &bank, int destination,
              int tag);

    // Receive.
    void receive(int num_recv, Fission_Site_Container_t &bank,
                 Fission_Site_Container_t &recv_bank, int destination,
                 profugus::Request &handle, int tag);

    // Number of sets and this set.
    int d_num_sets, d_set;

    // Left/right neighbors and the number of neighbors (at most 2) in set
    // space.
    int d_left, d_right, d_num_nbors;

    // Global number of fission sites sampled in this cycle (across all sets).
    int d_num_global;

    // The target number of fission sites on this process and set after
    // rebalance.
    int d_target_set;

    // Number of fission sites on each set.
    Vec_Int d_sites_set;

    // Current global fission bank array bounds on this set.
    Array_Bnds d_bnds;

    // Current global fission bank array bounds for the left and right
    // neighbors.
    Array_Bnds d_left_bnds, d_right_bnds;

    // Target global fission bank array bounds on this set.
    Array_Bnds d_target_bnds;

    // Receive banks.
    Fission_Site_Container_t d_recv_left, d_recv_right;

    // Counters for send/receive diagnostics.
    int d_num_recv, d_num_send, d_num_iter;

    // Receive handles.
    profugus::Request d_handle_left, d_handle_right;

    // Size of a fission site in bytes.
    int d_size_fs;
};

} // end namespace profugus

#endif // mc_Fission_Rebalance_hh

//---------------------------------------------------------------------------//
//              end of mc/Fission_Rebalance.hh
//---------------------------------------------------------------------------//
