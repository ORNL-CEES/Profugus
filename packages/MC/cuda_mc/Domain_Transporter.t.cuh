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
Domain_Transporter<Geometry>::Domain_Transporter(RCP_Std_DB     db,
                                                 SDP_Geometry   geometry,
                                                 SDP_Physics    physics,
                                                 SDP_VR         vr)
    : d_sample_fission_sites(false)
    , d_keff(0.0)
{
    REQUIRE(geometry.get_host_ptr());
    REQUIRE(geometry.get_device_ptr());
    REQUIRE(physics.get_host_ptr());
    REQUIRE(physics.get_device_ptr());
    d_geometry = geometry.get_device_ptr();
    d_physics  = physics.get_device_ptr();

    d_max_steps = db->get("max_steps",std::numeric_limits<int>::max());

    // Initialize tallier to null
    d_tallier = nullptr;

    // VR is optional
    if( vr.get_host_ptr() )
    {
        d_vr = vr.get_device_ptr();
        REQUIRE( d_vr );
    }
    else
    {
        d_vr = nullptr;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set tallier
 */
template <class Geometry>
void Domain_Transporter<Geometry>::set(SDP_Tallier tallier)
{
    REQUIRE( tallier.get_device_ptr() );
    d_tallier = tallier.get_device_ptr();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set fission site sampling.
 *
 * \param fission_sites
 * \param keff
 */
template <class Geometry>
void Domain_Transporter<Geometry>::set(SP_Fission_Site_Vec fission_sites,
                                       double              keff)
{
    // assign the container
    d_fission_sites = fission_sites->data().get();

    // Number of allocated sites
    d_max_fission_sites = fission_sites->size();

    // initialize the sampling flag
    d_sample_fission_sites = false;

    // set the flag indicating whether fission sites should be sampled or not
    if (d_max_fission_sites > 0)
    {
        d_sample_fission_sites = true;
    }

    // assign current iterate of keff
    d_keff = keff;

    // initialize the number of fission sites to 0
    d_num_fission_sites = 0;
}

} // end namespace cuda_mc

#endif // cuda_mc_Domain_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.t.cuh
//---------------------------------------------------------------------------//
