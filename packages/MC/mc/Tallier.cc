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

    d_geometry = geometry;
    d_physics  = physics;
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
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally any source events.
 */
void Tallier::source(const Particle_t &p)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform all end-history tally tasks.
 */
void Tallier::end_history()
{
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Tallier.cc
//---------------------------------------------------------------------------//
