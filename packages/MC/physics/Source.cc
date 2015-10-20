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

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR/DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Source::Source(SP_Geometry    geometry,
               SP_Physics     physics)
    : b_geometry(geometry)
    , b_physics(physics)
    , b_node(profugus::node())
    , b_nodes(profugus::nodes())
{
    REQUIRE(b_geometry);
    REQUIRE(b_physics);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Virtual destructor.
 */
Source::~Source()
{
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source.cc
//---------------------------------------------------------------------------//
