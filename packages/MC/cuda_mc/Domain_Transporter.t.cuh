//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter.t.cuh
 * \author Thomas M. Evans
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Domain_Transporter_t_cuh
#define cuda_mc_Domain_Transporter_t_cuh

#include <cmath>
#include <memory>

#include "harness/DBC.hh"
#include "harness/Diagnostics.hh"
#include "geometry/Definitions.hh"
#include "mc/Definitions.hh"
#include "Domain_Transporter.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Domain_Transporter<Geometry>::Domain_Transporter()
    : d_sample_fission_sites(false)
    , d_keff(0.0)
{
    d_roulette = false;
}


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the geometry and physics classes.
 *
 * \param geometry
 * \param physics
 */
template <class Geometry>
void Domain_Transporter<Geometry>::set(SDP_Geometry geometry,
                                       SDP_Physics  physics)
{
    REQUIRE(geometry.get_host_ptr());
    REQUIRE(geometry.get_device_ptr());
    REQUIRE(physics.get_host_ptr());
    REQUIRE(physics.get_device_ptr());

    // Get device pointers to Geometry and Physics
    d_geometry = geometry.get_device_ptr();
    d_physics  = physics.get_device_ptr();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the variance reduction
 *
 * \param vr
 */
template <class Geometry>
void Domain_Transporter<Geometry>::set(SDP_VR vr)
{
    REQUIRE( vr.get_host_ptr() );
    REQUIRE( vr.get_device_ptr() );

    // Get device pointer to VR
    d_roulette = true;

    d_vr = vr.get_device_ptr();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set regular tallies.
 *
 * \param tallies
 */
/*
template <class Geometry>
void Domain_Transporter<Geometry>::set(SP_Tallier tallies)
{
    REQUIRE(tallies);
    d_tallier = tallies;
    ENSURE(d_tallier);
}
*/

//---------------------------------------------------------------------------//
/*!
 * \brief Set fission site sampling.
 *
 * \param fission_sites
 * \param keff
 */
/*
template <class Geometry>
void Domain_Transporter<Geometry>::set(SP_Fission_Sites fission_sites,
                                       double           keff)
{
    // assign the container
    d_fission_sites = fission_sites;

    // initialize the sampling flag
    d_sample_fission_sites = false;

    // set the flag indicating whether fission sites should be sampled or not
    if (d_fission_sites)
    {
        d_sample_fission_sites = true;
    }

    // assign current iterate of keff
    d_keff = keff;

    // initialize the number of fission sites to 0
    d_num_fission_sites = 0;
}
*/

} // end namespace cuda_mc

#endif // cuda_mc_Domain_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.t.cuh
//---------------------------------------------------------------------------//
