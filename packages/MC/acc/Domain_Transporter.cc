//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Domain_Transporter.cc
 * \author Thomas M. Evans
 * \date   Sat Nov 01 10:56:35 2014
 * \brief  Domain_Transporter member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Domain_Transporter.hh"

namespace acc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Domain_Transporter::Domain_Transporter()
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS (called from the CPU)
//---------------------------------------------------------------------------//
/*!
 * \brief Set the geometry.
 *
 * This takes a CPU geometry and creates a device geometry from it.  The
 * device geometry is automatically allocated on the device.
 */
void Domain_Transporter::set(SP_Geometry geometry)
{
    // build the device geometry
    d_geometry = std::make_shared<Geometry>(*geometry);
}

} // end namespace acc

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.cc
//---------------------------------------------------------------------------//
