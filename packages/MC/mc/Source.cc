//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source.cc
 * \author Thomas M. Evans
 * \date   Mon May 05 14:28:41 2014
 * \brief  Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Source.hh"
#include "harness/DBC.hh"
#include "comm/global.hh"
#include "Global_RNG.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR/DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Source::Source(SP_Geometry    geometry,
               SP_Physics     physics,
               SP_RNG_Control rng_control)
    : b_geometry(geometry)
    , b_physics(physics)
    , b_rng_control(rng_control)
    , b_node(profugus::node())
    , b_nodes(profugus::nodes())
    , d_rng_stream(0)
{
    Require (b_geometry);
    Require (b_physics);
    Require (b_rng_control);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Virtual destructor.
 */
Source::~Source()
{
}

//---------------------------------------------------------------------------//
// PROTECTED FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate offsets for random numbers.
 */
void Source::make_RNG()
{
    // calculate offsets for generating random numbers on each processor (so
    // we don't use the same rng streams)

    // we use the same RNG for every particle in the cycle, each cycle a new
    // RNG is generated on each domain

    // make the random number generator on this domain for this cycle
    profugus::Global_RNG::d_rng = b_rng_control->rng(d_rng_stream + b_node);

    // advance to the next set of streams
    d_rng_stream += b_nodes;

    Ensure (profugus::Global_RNG::d_rng.assigned());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source.cc
//---------------------------------------------------------------------------//
