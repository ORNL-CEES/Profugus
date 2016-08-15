//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.t.hh
 * \author Stuart Slattery
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Tallier_t_cuh
#define cuda_mc_Tallier_t_cuh

#include <algorithm>
#include <utility>

#include "harness/Warnings.hh"
#include "cuda_utils/CudaDBC.hh"
#include "comm/Timing.hh"
#include "Tallier.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Tallier<Geometry>::Tallier()
    : d_build_phase(CONSTRUCTED)
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Add a pathlength tally.
 */
template <class Geometry>
void Tallier<Geometry>::add_pathlength_tally(SP_Pathlength_Tally tally)
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
template <class Geometry>
void Tallier<Geometry>::add_source_tally(SP_Source_Tally tally)
{
    REQUIRE(tally);
    REQUIRE(d_build_phase < BUILT);

    // add the tally
    d_src.push_back(tally);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize internal data structures after adding tallies.
 */
template <class Geometry>
void Tallier<Geometry>::build()
{
    REQUIRE(d_build_phase == CONSTRUCTED);
    REQUIRE(d_tallies.empty());

    // add pathlength and source tallies to the "totals"
    d_tallies.insert(d_tallies.end(), d_pl.begin(), d_pl.end());
    d_tallies.insert(d_tallies.end(), d_src.begin(), d_src.end());
    CHECK(num_tallies() == num_source_tallies() + num_pathlength_tallies() );

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
 */
template <class Geometry>
void Tallier<Geometry>::path_length(
    const cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    REQUIRE(d_build_phase == BUILT);

    if (!num_pathlength_tallies())
        return;

    SCOPED_TIMER_3("CUDA_MC::Tallier.path_length");

    // accumulate results for all pathlength tallies
    for (auto t : d_pl)
    {
        t->accumulate( particles );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally any source events.
 *
 * \param p particle
 */
template <class Geometry>
void Tallier<Geometry>::source(
    const cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    REQUIRE(d_build_phase == BUILT);

    if (!num_source_tallies())
        return;

    SCOPED_TIMER_3("CUDA_MC::Tallier.source");

    // accumulate results for all pathlength tallies
    for (auto t : d_src)
    {
        t->birth( particles );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin active kcode cycles.
 */
template <class Geometry>
void Tallier<Geometry>::begin_active_cycles()
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("CUDA_MC::Tallier.begin_active_cycles");

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
template <class Geometry>
void Tallier<Geometry>::begin_cycle()
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("CUDA_MC::Tallier.begin_cycle");

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
template <class Geometry>
void Tallier<Geometry>::end_cycle(double num_particles)
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("CUDA_MC::Tallier.end_cycle");

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
template <class Geometry>
void Tallier<Geometry>::end_history()
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("CUDA_MC::Tallier.end_history");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->end_history();
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
template <class Geometry>
void Tallier<Geometry>::finalize(double num_particles)
{
    REQUIRE(d_build_phase == BUILT);

    SCOPED_TIMER_2("CUDA_MC::Tallier.finalize");

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
template <class Geometry>
void Tallier<Geometry>::reset()
{
    REQUIRE(d_build_phase == FINALIZED);

    SCOPED_TIMER_2("CUDA_MC::Tallier.reset");

    // begin active for each tally
    for (auto t : d_tallies)
    {
        t->reset();
    }

    // clear the list of tallies (need to call build again to get these)
    d_tallies.clear();

    // set the build phase
    d_build_phase = CONSTRUCTED;

    ENSURE(!is_finalized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two talliers.
 *
 * This is useful for temporarily deactivating tallying (say, during inactive
 * cycles in a kcode calculation).
 */
template <class Geometry>
void Tallier<Geometry>::swap(Tallier<Geometry> &rhs)
{
    // swap vector internals
    d_pl.swap(rhs.d_pl);
    d_src.swap(rhs.d_src);
    d_tallies.swap(rhs.d_tallies);

    // swap build phase
    std::swap(d_build_phase, rhs.d_build_phase);
}

} // end namespace cuda_profugus

#endif // cuda_mc_Tallier_t_cuh

//---------------------------------------------------------------------------//
//                 end of Tallier.t.hh
//---------------------------------------------------------------------------//
