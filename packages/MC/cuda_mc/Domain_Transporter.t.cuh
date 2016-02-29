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
#include "utils/String_Functions.hh"
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
Domain_Transporter<Geometry>::Domain_Transporter(RCP_Std_DB db)
    : d_sample_fission_sites(false)
    , d_keff(0.0)
{
    // Build variance reduction
    std::string var = profugus::to_lower(
        db->get<std::string>("variance reduction",std::string("roulette")) );
    if( var == "roulette" )
    {
        std::cout << "Building DT with Roulette" << std::endl;
        d_roulette = true;
        auto sp_vr = std::make_shared<VR_Roulette_t>(db);
        d_vr_host = SDP_VR( sp_vr );
        d_vr = d_vr_host.get_device_ptr();
    }
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
    d_geometry_host = geometry;
    d_geometry = d_geometry_host.get_device_ptr();

    d_physics_host  = physics;
    d_physics = d_physics_host.get_device_ptr();
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
