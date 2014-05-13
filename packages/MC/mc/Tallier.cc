//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tallier.cc
 * \author Thomas M. Evans
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "Tallier.hh"

namespace profugus
{

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
    Require (geometry);
    Require (physics);
    Require (d_build_phase == CONSTRUCTED);

    d_geometry = geometry;
    d_physics  = physics;

    // set the build phase
    d_build_phase = ASSIGNED;

    Ensure (d_build_phase == ASSIGNED);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize internal data structures after adding tallies.
 */
void Tallier::build()
{
    Require (d_build_phase == ASSIGNED);

    // Set the build phase
    d_build_phase = BUILT;

    Ensure (d_build_phase == BUILT);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process path-length tally events.
 *
 * \param step
 * \param p
 */
void Tallier::path_length(double            step,
                          const Particle_t &p)
{
    Require (d_build_phase == BUILT);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally any source events.
 */
void Tallier::source(const Particle_t &p)
{
    Require (d_build_phase == BUILT);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform all end-history tally tasks.
 */
void Tallier::end_history()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize tallies.
 */
void Tallier::finalize(int Np)
{
    Require (d_build_phase == BUILT);

    // set the build phase
    d_build_phase = FINALIZED;

    Ensure (is_finalized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset results in all tallies
 *
 * This does not remove tallies from the tallier: it will just reset the
 * accumulators. See implementation of the Tally daughter classes for details,
 * but generally this doesn't clear the values in existing smart-pointer fields.
 *
 * \pre \c finalize() was called on tallies
 * \post \c is_finalized() returns false
 */
void Tallier::reset()
{
    Require (d_build_phase == FINALIZED);

    // set the build phase
    d_build_phase = ASSIGNED;

    Ensure (!is_finalized());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Tallier.cc
//---------------------------------------------------------------------------//
