//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tallier.cc
 * \author Thomas M. Evans
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <utility>

#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "comm/Timing.hh"
#include "Tallier.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// PRIVATE TEMPLATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Prune added tallies for doubles.
 */
template<class Vec_T>
void Tallier::prune(Vec_T &tallies)
{
    // sort the tally container based on the tally name
    auto sort_f = [](const SP_Tally &a, const SP_Tally &b)
                  { return a->name() < b->name(); };
    std::sort(tallies.begin(), tallies.end(), sort_f);

    // make a new tally container
    Vec_T new_tallies;

    // define a place-holder for the "previous tally" to use to check for
    // duplicates
    typename Vec_T::value_type prev_tally;
    CHECK(!prev_tally);

    // iterate through tallies and remove duplicates
    for (const auto &t : tallies)
    {
        CHECK(t);

        // if this is a dublicate, continue
        if (t == prev_tally)
            continue;

        // add the tally
        new_tallies.push_back(t);

        // update the previous tally
        prev_tally = t;
    }
    CHECK(new_tallies.size() <= tallies.size());
    REMEMBER(int size = new_tallies.size());

    // swap the new_tallies with the tallies
    tallies.swap(new_tallies);
    ENSURE(tallies.size() == size);
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Tallier::Tallier()
    : d_build_phase(CONSTRUCTED)
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the geometry and physics.
 *
 * \param geometry
 * \param physics
 */
void Tallier::set(SP_Geometry geometry,
                  SP_Physics  physics)
{
    REQUIRE(geometry);
    REQUIRE(physics);
    REQUIRE(d_build_phase == CONSTRUCTED);

    d_geometry = geometry;
    d_physics  = physics;

    // set the build phase
    d_build_phase = ASSIGNED;

    ENSURE(d_build_phase == ASSIGNED);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a pathlength tally.
 */
void Tallier::add_pathlength_tally(SP_Pathlength_Tally tally)
{
    REQUIRE(tally);
    REQUIRE(d_build_phase < BUILT);

    // add the tally
    d_pl.push_back(tally);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a source tally.
 */
void Tallier::add_source_tally(SP_Source_Tally tally)
{
    REQUIRE(tally);
    REQUIRE(d_build_phase < BUILT);

    // add the tally
    d_src.push_back(tally);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a compound tally.
 */
void Tallier::add_compound_tally(SP_Compound_Tally tally)
{
    REQUIRE(tally);
    REQUIRE(d_build_phase < BUILT);

    // add the compound tally
    d_comp.push_back(tally);

    // add the source and pathlength tallies owned by the compound tally
    add_source_tally(tally->get_src_tally());
    add_pathlength_tally(tally->get_pl_tally());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize internal data structures after adding tallies.
 */
void Tallier::build()
{
    REQUIRE(d_build_phase == ASSIGNED);
    REQUIRE(d_tallies.empty());

    // prune the tallies for duplicates
    prune(d_pl);
    prune(d_src);
    prune(d_comp);

    // add pathlength, source, and compound tallies to the "totals"
    d_tallies.insert(d_tallies.end(), d_pl.begin(), d_pl.end());
    d_tallies.insert(d_tallies.end(), d_src.begin(), d_src.end());
    d_tallies.insert(d_tallies.end(), d_comp.begin(), d_comp.end());
    CHECK(num_tallies() == num_source_tallies() + num_pathlength_tallies() +
          num_compound_tallies());

    // Set the build phase
    d_build_phase = BUILT;

    // Print warnings if applicable
    if (num_tallies() == 0)
    {
        ADD_WARNING("No tallies are being used.");
    }

    ENSURE(d_build_phase == BUILT);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process path-length tally events.
 *
 * \param step step-length
 * \param p particle
 */
void Tallier::path_length(double            step,
                          const Particle_t &p)
{
    REQUIRE(d_build_phase == BUILT);
    REQUIRE(step >= 0.0);

    if (!num_pathlength_tallies())
        return;

    SCOPED_TIMER_3("MC::Tallier.path_length");

    // accumulate results for all pathlength tallies
    for (auto t : d_pl)
    {
        t->accumulate(step, p);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally any source events.
 *
 * \param p particle
 */
void Tallier::source(const Particle_t &p)
{
    REQUIRE(d_build_phase == BUILT);

    if (!num_source_tallies())
        return;

    SCOPED_TIMER_3("MC::Tallier.source");

    // accumulate results for all pathlength tallies
    for (auto t : d_src)
    {
        t->birth(p);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin active kcode cycles.
 */
void Tallier::begin_active_cycles()
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("MC::Tallier.begin_active_cycles");

    // begin active cycles for each tally
    for (auto t : d_tallies)
    {
        t->begin_active_cycles();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin a new cycle in a kcode calculation.
 */
void Tallier::begin_cycle()
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("MC::Tallier.begin_cycle");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->begin_cycle();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to end a cycle in a kcode calculation.
 */
void Tallier::end_cycle(double num_particles)
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("MC::Tallier.end_cycle");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->end_cycle(num_particles);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform all end-history tally tasks.
 */
void Tallier::end_history(const Particle_t &p)
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("MC::Tallier.end_history");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->end_history(p);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize tallies.
 *
 * Does post-solve processing of tallies including parallel reductions and
 * normalization.
 *
 * \post is_finalized() == true
 */
void Tallier::finalize(double num_particles)
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("MC::Tallier.finalize");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->finalize(num_particles);
    }

    // set the build phase
    d_build_phase = FINALIZED;

    ENSURE(is_finalized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset results in all tallies
 *
 * This does not remove tallies from the tallier: it will just reset the
 * accumulators. See implementation of the Tally daughter classes for details,
 * but generally this doesn't clear the values in existing smart-pointer
 * fields.
 *
 * \pre \c finalize() was called on tallies
 * \post \c is_finalized() returns false
 */
void Tallier::reset()
{
    REQUIRE(d_build_phase == FINALIZED);

    SCOPED_TIMER_2("MC::Tallier.reset");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->reset();
    }

    // clear the list of tallies (need to call build again to get these)
    d_tallies.clear();

    // set the build phase
    d_build_phase = ASSIGNED;

    ENSURE(!is_finalized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two talliers.
 *
 * This is useful for temporarily deactivating tallying (say, during inactive
 * cycles in a kcode calculation).
 */
void Tallier::swap(Tallier &rhs)
{
    // swap vector internals
    d_pl.swap(rhs.d_pl);
    d_src.swap(rhs.d_src);
    d_comp.swap(rhs.d_comp);
    d_tallies.swap(rhs.d_tallies);

    // swap geometry and physics
    std::swap(d_geometry, rhs.d_geometry);
    std::swap(d_physics, rhs.d_physics);

    // swap build phase
    std::swap(d_build_phase, rhs.d_build_phase);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Tallier.cc
//---------------------------------------------------------------------------//
