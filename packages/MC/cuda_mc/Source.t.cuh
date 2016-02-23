//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source.t.cuh
 * \author Steven Hamilton
 * \date   Mon May 05 14:28:41 2014
 * \brief  Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_t_cuh
#define cuda_mc_Source_t_cuh

#include "Source.cuh"
#include "harness/DBC.hh"
#include "comm/global.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR/DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Source<Geometry>::Source(SDP_Geometry geometry)
    : b_node(profugus::node())
    , b_nodes(profugus::nodes())
    , b_geometry_host(geometry)
    , d_rng_stream(0)
{
    b_geometry = b_geometry_host.get_device_ptr();

    REQUIRE(b_geometry);
}

//---------------------------------------------------------------------------//
// PROTECTED FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate offsets for random numbers.
 */
template <class Geometry>
void Source<Geometry>::make_RNG()
{
    // calculate offsets for generating random numbers on each processor (so
    // we don't use the same rng streams)

    // we use the same RNG for every particle in the cycle, each cycle a new
    // RNG is generated on each domain

    // TODO: What do we want to do here to set up the on-device RNG?

    // make the random number generator on this domain for this cycle
    //profugus::Global_RNG::d_rng = b_rng_control->rng(d_rng_stream + b_node);

    // advance to the next set of streams
    //d_rng_stream += b_nodes;

    //ENSURE(profugus::Global_RNG::d_rng.assigned());
}

} // end namespace profugus

#endif // cuda_mc_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source.t.cuh
//---------------------------------------------------------------------------//
