//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Definitions.hh
 * \author Thomas M. Evans
 * \date   Mon Aug 22 11:21:08 2011
 * \brief  General geometry functions and types.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_Definitions_hh
#define geometry_Definitions_hh

#include <cmath>
#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "utils/Definitions.hh"

namespace denovo
{

namespace geometry
{

typedef unsigned int matid_type;
typedef unsigned int cell_type;

//---------------------------------------------------------------------------//
// ENUMERATED TYPES
//---------------------------------------------------------------------------//
/*!
 * \enum Boundary_State
 * \brief Current state with respect to outer geometry boundary
 *
 */
enum Boundary_State {
    INSIDE  = 0, //!< Particle is within geometry
    OUTSIDE,     //!< Particle is outside geometry
    REFLECT      //!< Particle on boundary surface and has just reflected
};

} // end namespace geometry

} // end namespace denovo

#endif // geometry_Definitions_hh

//---------------------------------------------------------------------------//
//              end of geometry/Definitions.hh
//---------------------------------------------------------------------------//
