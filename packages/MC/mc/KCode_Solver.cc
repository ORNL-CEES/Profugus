//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/KCode_Solver.cc
 * \author Thomas M. Evans
 * \date   Mon May 19 10:30:32 2014
 * \brief  KCode_Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KCode_Solver.hh"

#include <iostream>
#include <iomanip>

#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
KCode_Solver::KCode_Solver(RCP_Std_DB db)
    : d_db(db)
    , d_build_phase(CONSTRUCTED)
    , d_quiet(db->get("quiet", false))
{
    VALIDATE(profugus::num_available_threads() == 1,
             "Cannot run multithreaded-eigenvalue problems.");

    // set quiet off on work nodes
    if (profugus::node() != 0)
        d_quiet = true;

    ENSURE(d_build_phase == CONSTRUCTED);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the underlying fixed-source transporter and fission source.
 */
void KCode_Solver::set(SP_Source_Transporter transporter,
                       SP_Fission_Source     source)
{
    REQUIRE(d_build_phase == CONSTRUCTED);
    REQUIRE(transporter);
    REQUIRE(source);

    // assign the solver
    d_transporter = transporter;

    // assign the source
    d_source = source;

    // fission site container
    d_fission_sites = d_source->create_fission_site_container();
    CHECK(d_fission_sites);
    CHECK(d_fission_sites->empty());

    // store the integrated cycle weight (const number of requested particles
    // per cycle)
    d_Np = static_cast<double>(d_source->Np());
    CHECK(d_Np > 0.0);

    // get a reference to the tallier
    b_tallier = d_transporter->tallier();
    INSIST(b_tallier, "The tallier has not been assigned.");
    CHECK(b_tallier->geometry() && b_tallier->physics());

    // get initial k and build keff tally
    double init_keff = d_db->get("keff_init", 1.0);
    VALIDATE(init_keff >= 0., "Initial keff guess (keff_init="
            << init_keff << ") must be nonnegative.");
    b_keff_tally = std::make_shared<Keff_Tally>(init_keff,
                                                b_tallier->physics());

    // create our "disabled" (inactive cycle) tallier
    d_inactive_tallier = std::make_shared<Tallier_t>();
    d_inactive_tallier->set(b_tallier->geometry(), b_tallier->physics());

    d_build_phase = ASSIGNED;

    ENSURE(b_tallier);
    ENSURE(b_keff_tally);
    ENSURE(d_transporter);
    ENSURE(d_build_phase == ASSIGNED);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the acceleration method.
 */
void KCode_Solver::set(SP_FM_Acceleration acceleration)
{
    d_acceleration = acceleration;
}

//---------------------------------------------------------------------------//
// INHERITED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Solve the kcode problem.
 *
 * This runs the initialize/iterate/finalize process, and also prints user
 * output.
 */
void KCode_Solver::solve()
{
    using std::endl; using std::cout;
    using std::setw; using std::fixed; using std::scientific;

    SCOPED_TIMER("MC::KCode_Solver.solve");

    // set cycle numbers
    const int num_total    = d_db->get("num_cycles",          50);
    const int num_inactive = d_db->get("num_inactive_cycles", 10);
    const int num_active   = num_total - num_inactive;

    VALIDATE(num_inactive > 0,
            "The number of  inactive keff cycles (num_inactive_cycles="
            << num_inactive << ") must be positive");
    VALIDATE(num_total > num_inactive,
            "The number of keff cycles (num_cycles=" << num_total << ") "
            "must be greater than the number of inactive cycles ("
            "num_inactive_cycles=" << num_inactive << ")");

    ENSURE(num_inactive + num_active == num_total);
    ENSURE(num_total > 0);

    // Preallocate diagnostics
    DIAGNOSTICS_ONE(vec_integers["np_fission"]      .reserve(num_total );)
    DIAGNOSTICS_TWO(vec_doubles ["local_np_fission"].reserve(num_total );)

    // Timers (some required, some optional)
    Timer cycle_timer;

    // Set up initial source etc.
    this->initialize();

    // Print column header
    if (!d_quiet)
    {
        cout << ">>> Beginning inactive cycles." << endl;
        cout << endl;
        cout << " Cycle  k_cycle   Time (s) " << endl;
    }

    // Solve inactive cycles
    for (int cycle = 0; cycle < num_inactive; ++cycle)
    {
        // Iterate and time
        cycle_timer.start();
        this->iterate();
        cycle_timer.stop();

        // Print diagnostic output
        if (!d_quiet)
        {
            cout.precision(6);
            cout.setf(std::ios::internal);

            cout << fixed      << setw(4) << cycle << "   "
                 << fixed      << setw(8) << b_keff_tally->latest() << "  "
                 << scientific << setw(9) << cycle_timer.TIMER_CLOCK()
                 << endl;
        }
    }

    // Prepare for active cycles
    CHECK(num_cycles() == num_inactive);
    this->begin_active_cycles();
    CHECK(num_cycles() == 0);

    // Print column header
    if (!d_quiet)
    {
        cout << endl;
        cout << ">>> Beginning active cycles." << endl;
        cout << endl;
        cout << " Cycle  k_cycle    k_avg    k_variance    Time (s) "
             << endl;
    }

    // Solve active cycles
    for (int cycle = 0; cycle < num_active; ++cycle)
    {
        // Iterate and time
        cycle_timer.start();
        this->iterate();
        cycle_timer.stop();

        // Print diagnostic output
        if (!d_quiet)
        {
            cout.precision(6);
            cout.setf(std::ios::internal);

            cout << fixed      << setw(4)  << cycle << "   "
                 << fixed      << setw(8)  << b_keff_tally->latest() << "  "
                 << fixed      << setw(8)  << b_keff_tally->mean() << "  "
                 << scientific << setw(11) << b_keff_tally->variance() << "  "
                 << scientific << setw(11) << cycle_timer.TIMER_CLOCK()
                 << endl;
        }
    }
    CHECK(num_cycles() == num_active);

    if (!d_quiet)
        cout << endl;

    profugus::global_barrier();

    // Finalize
    this->finalize();

    ENSURE(b_tallier->is_finalized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Call to reset the solver and tallies for another kcode run
 */
void KCode_Solver::reset()
{
    INSIST(d_build_phase == FINALIZED,
            "Reset may be called only after finalizing.");

    // Reset tallies
    b_tallier->reset();
    d_inactive_tallier->reset();

    d_build_phase = ASSIGNED;

    ENSURE(!b_tallier->is_built());
    ENSURE(d_build_phase == ASSIGNED);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the k-eigenvalue problem.
 *
 * The initialize method adds the keffective tally to the active/inactive
 * talliers, builds the initial/fission source, and changes the tallier being
 * used.
 *
 * This last step allows user-added tallies (e.g. flux tallies) to be turned
 * off during the initial cycles, so that their biased distributions don't
 * contribute to the final (active-averaged) results.
 *
 * After this "swap", the keff tally (and potentially other diagnostic tallies
 * added to the 'inactive' tallier) are the only ones that will be called by the
 * transporter.
 */
void KCode_Solver::initialize()
{
    INSIST(d_build_phase == ASSIGNED,
            "initialize must be called only after calling set()");

    INSIST(!b_tallier->is_built(), "The given tallier can only be built once.");
    INSIST(!b_tallier->is_finalized(), "The given tallier was used in "
            "a prior transport solve without reset() being called.");

    // Add the keff tally to the ACTIVE cycle tallier and build it
    b_tallier->add_pathlength_tally(b_keff_tally);
    b_tallier->build();
    CHECK(b_keff_tally->inactive_cycle_tally());
    CHECK(b_tallier->is_built());
    CHECK(b_tallier->num_pathlength_tallies() >= 1);

    // get tallies from the active tallier and add them to the inactive
    // tallier if they are inactive tallies
    for (auto titr = b_tallier->begin(); titr != b_tallier->end(); ++titr)
    {
        CHECK(*titr);

        // get the SP to the tally
        auto tally = *titr;

        // if this tally should be on during inactive cycles, add it
        if (tally->inactive_cycle_tally())
        {
            // create tallies (only 1 can be valid)
            auto pl_t  = std::dynamic_pointer_cast<Pathlength_Tally>(tally);
            auto src_t = std::dynamic_pointer_cast<Source_Tally>(tally);
            auto cpd_t = std::dynamic_pointer_cast<Compound_Tally>(tally);

            // attempt to cast to valid tally types and add the tally
            if (pl_t)
            {
                CHECK(!src_t && !cpd_t);
                d_inactive_tallier->add_pathlength_tally(pl_t);
            }
            else if (src_t)
            {
                CHECK(!pl_t && !cpd_t);
                d_inactive_tallier->add_source_tally(src_t);
            }
            else if (cpd_t)
            {
                CHECK(!pl_t && !src_t);
                d_inactive_tallier->add_compound_tally(cpd_t);
            }
            else
            {
                throw profugus::assertion("Unknown tally type.");
            }
        }
    }

    // build the inactive tallies
    d_inactive_tallier->build();
    CHECK(d_inactive_tallier->is_built());
    CHECK(d_inactive_tallier->num_pathlength_tallies() >= 1);

    // Swap our temp tallier with the user-specified tallies so that we
    // initially only tally k effective.
    // b_tallier is always the one actively getting called by transporter
    swap(*d_inactive_tallier, *b_tallier);

    // initialize the acceleration
    if (d_acceleration)
    {
        d_acceleration->initialize(d_db);

        // build the initial fission source (which may or may not use an
        // initial distribution depending on the acceleration options)
        if (d_source->is_initial_source())
        {
            d_acceleration->build_initial_source(*d_source);
        }
    }
    // build the initial source if no acceleration is on
    else
    {
        // build the initial fission source
        if (d_source->is_initial_source())
        {
            d_source->build_initial_source();
        }
    }

    d_build_phase = INACTIVE_SOLVE;

    ENSURE(b_tallier->num_pathlength_tallies() >= 1);
    ENSURE(b_keff_tally->cycle_count() == 0);
    ENSURE(d_build_phase == INACTIVE_SOLVE);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform one k-eigenvalue cycle.
 */
void KCode_Solver::iterate()
{
    REQUIRE(d_fission_sites && d_fission_sites->empty());
    REQUIRE(b_tallier->is_built() && !b_tallier->is_finalized());
    REQUIRE(d_build_phase == INACTIVE_SOLVE || d_build_phase == ACTIVE_SOLVE);

    // store the number of source particles for this cycle
    DIAGNOSTICS_ONE(vec_integers["np_fission"].push_back(
                        d_source->total_num_to_transport()));
    DIAGNOSTICS_TWO(vec_integers["local_np_fission"].push_back(
                        d_source->num_to_transport()));

    // assign the current source state in the transporter (this will generally
    // be a pass through, but it gives the solver a chance to update
    // quantities if the source changes from cycle-to-cycle)
    d_transporter->assign_source(d_source);

    // set the solver to sample fission sites
    d_transporter->sample_fission_sites(d_fission_sites,
                                        b_keff_tally->latest());

    // inititalize the acceleration at the beginning of the cycle with the
    // rebalanced fission sites owned by the source
    if (d_acceleration)
    {
        d_acceleration->start_cycle(b_keff_tally->latest(),
                                    d_source->fission_sites());
    }

    // initialize keff tally to the beginning of the cycle
    b_tallier->begin_cycle();

    // solve the fixed source problem using the transporter
    d_transporter->solve();

    // do end-of-cycle tally processing including global sum Note: this is the
    // total *requested* number of particles, not the actual number of
    // histories. Each processor adjusts the particle weights so that the
    // total weight emitted, summed over all processors, is d_Np.
    b_tallier->end_cycle(d_Np);

    // process acceleration at the end of the cycle
    if (d_acceleration)
    {
        d_acceleration->end_cycle(*d_fission_sites);

        // update the number of fission sites (source particles) per cycle
        if (d_acceleration->update_num_particles())
        {
            d_Np = d_acceleration->num_sites();

            // update the source
            d_source->update_Np(d_Np);
        }
    }

    // build a new source from the fission site distribution
    d_source->build_source(d_fission_sites);

    ENSURE(d_fission_sites);
    ENSURE(d_fission_sites->empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Call to transition from inactive to active cycles
 *
 * This reinstates the original tallies and calls "finalize" on the inactive
 * tallies.
 */
void KCode_Solver::begin_active_cycles()
{
    REQUIRE(b_tallier->is_built() && !b_tallier->is_finalized());
    REQUIRE(d_inactive_tallier);

    INSIST(d_build_phase == INACTIVE_SOLVE,
            "begin_active_cycles must be called only after "
            "initializing and iterating on inactive cycles.");

    // Swap the saved user-specified tallies with the inactive-cycle tallier
    // so that we start tallying all the other functions
    swap(*d_inactive_tallier, *b_tallier);

    // Finalize the inactive tallier (will set the state to FINALIZED, but
    // shouldn't do anything else unless the user has explicitly added a tally
    // to it. This ensures that if the user has added, say, the same mesh
    // tally to both inactive and active cycles, they get an obvious "tally
    // was already finalized" error instead of a more subtle error (viz.,
    // being normalized to the *active* particle count instead of the *active
    // + inactive* particle count).
    d_inactive_tallier->finalize(num_cycles() * d_Np);

    // Tell all tallies to begin tallying active cycles
    b_tallier->begin_active_cycles();

    d_build_phase = ACTIVE_SOLVE;

    ENSURE(d_inactive_tallier->is_finalized());
    ENSURE(num_cycles() == 0);
    ENSURE(d_build_phase == ACTIVE_SOLVE);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize the k-eigenvalue problem.
 *
 * This calls "finalize" on the active tallies, doing global sums,
 * redistributions, etc.
 */
void KCode_Solver::finalize()
{
    INSIST(d_build_phase == ACTIVE_SOLVE,
            "Finalize can only be called after iterating with active cycles.");
    INSIST(num_cycles() > 0, "No active cycles were performed.");

    // Finalize tallies using global number particles
    CHECK(!b_tallier->is_finalized());
    b_tallier->finalize(num_cycles() * d_Np);

    d_build_phase = FINALIZED;

    ENSURE(b_tallier->is_finalized());
    ENSURE(d_build_phase == FINALIZED);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of KCode_Solver.cc
//---------------------------------------------------------------------------//
