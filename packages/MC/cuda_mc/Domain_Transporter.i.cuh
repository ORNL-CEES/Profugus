//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter.i.cuh
 * \author Thomas M. Evans
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Domain_Transporter_i_hh
#define cuda_mc_Domain_Transporter_i_hh

#include "cuda_utils/CudaDBC.hh"
#include "geometry/Definitions.hh"
#include "mc/Definitions.hh"
#include "Domain_Transporter.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Transport a particle through the domain.
 *
 * \param particle
 * \param bank
 */
template <class Geometry>
__device__
void Domain_Transporter<Geometry>::transport(Particle_t &particle) const
{
    REQUIRE(d_geometry);
    REQUIRE(d_physics);
    REQUIRE(particle.alive());

    // particle state
    Geo_State_t &geo_state = particle.geo_state();

    Step_Selector step;
    double dist_mfp, dist_col, dist_bnd, xs_tot;

    // step particle through domain while particle is alive; life is relative
    // to the domain, so a particle leaving the domain would be no longer
    // alive wrt the current domain even though the particle might be alive to
    // another domain
    while (particle.alive())
    {
        // calculate distance to collision in mean-free-paths
        dist_mfp = -log(particle.ran());

        // while we are hitting boundaries, continue to transport until we get
        // to the collision site
        particle.set_event(profugus::events::BOUNDARY);

        // process particles through internal boundaries until it hits a
        // collision site or leaves the domain
        //
        // >>>>>>> IMPORTANT <<<<<<<<
        // it is absolutely essential to check the geometry distance first,
        // before the distance to boundary mesh; this naturally takes care of
        // coincident boundary mesh and problem geometry surfaces (you may end
        // up with a small, negative distance to boundary mesh on the
        // subsequent step, but who cares)
        // >>>>>>>>>>>>-<<<<<<<<<<<<<
        while (particle.event() == profugus::events::BOUNDARY)
        {
            CHECK(particle.alive());
            CHECK(dist_mfp > 0.0);

            // total interaction cross section
            xs_tot = d_physics->total(profugus::physics::TOTAL, particle);
            CHECK(xs_tot >= 0.0);

            // sample distance to next collision
            if (xs_tot > 0.0)
                dist_col = dist_mfp / xs_tot;
            else
                dist_col = profugus::constants::huge;

            // initialize the distance to collision in the step selector
            step.initialize(dist_col, profugus::events::COLLISION);

            // calculate distance to next geometry boundary
            dist_bnd = d_geometry->distance_to_boundary(geo_state);
            step.submit(dist_bnd, profugus::events::BOUNDARY);

            // set the next event in the particle
            CHECK(step.tag() < profugus::events::END_EVENT);
            particle.set_event(static_cast<profugus::events::Event>(
                step.tag()));

            // path-length tallies (the actual movement of the particle will
            // take place when we process the various events)
            if( d_tallier )
                d_tallier->path_length(step.step(), particle);

            // update the mfp distance travelled
            dist_mfp -= step.step() * xs_tot;

            // process a particle through the geometric boundary
            if (particle.event() == profugus::events::BOUNDARY)
                process_boundary(particle);
        }

        // we are done moving to a problem boundary, now process the
        // collision; the particle is actually moved from its current position
        // in each of these calls

        // process particle at a collision
        if (particle.event() == profugus::events::COLLISION)
            process_collision(particle,step.step());

        // any future events go here ...
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

template <class Geometry>
__device__ void
Domain_Transporter<Geometry>::process_boundary(Particle_t &particle) const
{
    REQUIRE(particle.alive());
    REQUIRE(particle.event() == profugus::events::BOUNDARY);

    // return if not a boundary event

    // move a particle to the surface and process the event through the surface
    d_geometry->move_to_surface(particle.geo_state());

    // get the in/out state of the particle
    int state = d_geometry->boundary_state(particle.geo_state());

    // reflected flag
    bool reflected = false;

    // process the boundary crossing based on the geometric state of the
    // particle
    switch (state)
    {
        case profugus::geometry::OUTSIDE:
            // the particle has left the problem geometry
            particle.set_event(profugus::events::ESCAPE);
            particle.kill();

            ENSURE(particle.event() == profugus::events::ESCAPE);
            break;

        case profugus::geometry::REFLECT:
            // the particle has hit a reflecting surface
            reflected = d_geometry->reflect(particle.geo_state());
            CHECK(reflected);

            ENSURE(particle.event() == profugus::events::BOUNDARY);
            break;

        case profugus::geometry::INSIDE:
            // otherwise the particle is at an internal geometry boundary;
            // update the material id of the region the particle has entered
            particle.set_matid(d_geometry->matid(particle.geo_state()));

            // add variance reduction at surface crossings
            if( d_vr )
                d_vr->post_surface(particle);

            ENSURE(particle.event() == profugus::events::BOUNDARY);
            break;

        default:
            CHECK(0);
    }
}

//---------------------------------------------------------------------------//

template <class Geometry>
__device__ void
Domain_Transporter<Geometry>::process_collision(Particle_t &particle,
                                                double      dist) const
{
    REQUIRE(particle.event() == profugus::events::COLLISION);

    // move the particle to the collision site
    d_geometry->move_to_point(dist, particle.geo_state());

    // sample fission sites
    if (d_sample_fission_sites)
    {
        CHECK(d_fission_sites);
        CHECK(d_keff > 0.0);
        int num_sites = d_physics->sample_fission_site(particle, d_keff);

        if( num_sites > 0 )
        {
            // Create fission site
            Fission_Site site;
            site.m = particle.matid();
            site.r = d_geometry->position(particle.geo_state());

            // Get index of first site to be written using an atomic
            int offset = atomicAdd(&d_num_fission_sites,num_sites);
            
            // Don't write past end of allocated array
            num_sites = min(num_sites,d_max_fission_sites-offset);

            // Add fission sites
            for( int i = 0; i < num_sites; ++i )
                d_fission_sites[offset+i] = site;
        }
    }

    // use the physics package to process the collision
    d_physics->collide(particle);

    // apply weight windows
    if( d_vr )
        d_vr->post_collision(particle);
}

} // end namespace cuda_mc

#endif // cuda_mc_Domain_Transporter_i_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.i.cuh
//---------------------------------------------------------------------------//
