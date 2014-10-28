//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/Fixed_Source_Solver.cc
 * \author Thomas M. Evans
 * \date   Tue May 13 14:40:06 2014
 * \brief  Fixed_Source_Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fixed_Source_Solver.hh"

#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/P_Stream.hh"
#include "comm/Timing.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Fixed_Source_Solver::Fixed_Source_Solver()
    : d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the underlying source transporter and source.
 */
void Fixed_Source_Solver::set(SP_Source_Transporter transporter,
                              SP_Source             source)
{
    REQUIRE(transporter);
    REQUIRE(source);

    // assign the solver
    d_transporter = transporter;

    // assign the source
    d_source = source;

    // set the global, total number of particles to run
    d_Np = d_source->total_num_to_transport();
    DIAGNOSTICS_ONE(vec_integers["np"].push_back(d_Np));

    // get the tallies and assign them
    b_tallier = d_transporter->tallier();
    INSIST(b_tallier,
            "Tally not assigned in Source_Transporter in fixed-source solver.");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 *
 * This also finalizes tallies.
 */
void Fixed_Source_Solver::solve()
{
    using profugus::endl; using profugus::pcout;

    REQUIRE(d_source);
    REQUIRE(b_tallier);
    INSIST(!b_tallier->is_built(), "The given tallier can only be built once.");
    INSIST(!b_tallier->is_finalized(), "The given tallier was used in "
           "a prior transport solve without reset() being called.");

    // Build the tallier
    b_tallier->build();
    CHECK(b_tallier->is_built());

    SCOPED_TIMER("MC::Fixed_Source_Solver.solve");

    // store the number of source particles
    DIAGNOSTICS_TWO(vec_integers["local_np"].push_back(
                        d_source->num_to_transport()));

    // assign the current source state in the source transporter
    d_transporter->assign_source(d_source);

    // start the timer
    profugus::Timer fixed_source_timer;
    fixed_source_timer.start();

    // solve the fixed source problem using the transporter
    d_transporter->solve();

    // stop the timer
    profugus::global_barrier();
    fixed_source_timer.stop();

    // output
    pcout << ">>> Finished transporting " << profugus::setw(8)
          << d_source->total_num_to_transport()  << " particles in "
          << profugus::fixed << profugus::setw(12)
          << profugus::setprecision(6) << fixed_source_timer.TIMER_CLOCK()
          << " seconds" << endl;

    // Finalize tallies using global number of particles
    b_tallier->finalize(d_Np);
    ENSURE(b_tallier->is_finalized());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.cc
//---------------------------------------------------------------------------//
