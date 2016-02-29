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
    : b_geometry( geometry.get_device_ptr() )
{
    REQUIRE(b_geometry);
}

} // end namespace profugus

#endif // cuda_mc_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source.t.cuh
//---------------------------------------------------------------------------//
