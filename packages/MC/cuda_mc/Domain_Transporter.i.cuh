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
void Domain_Transporter<Geometry>::transport(int                pid,
                                             Particle_Vector_t &particles) const
{
    DEVICE_REQUIRE(d_geometry);
    DEVICE_REQUIRE(d_physics);
    DEVICE_REQUIRE(particles.alive(pid));

    // particle state
    Geo_State_t &geo_state = particles.geo_state(pid);

    Step_Selector step;
    double dist_mfp, dist_col, dist_bnd, xs_tot;

    int num_steps = 0;

    // step particle through domain while particle is alive; life is relative
    // to the domain, so a particle leaving the domain would be no longer
    // alive wrt the current domain even though the particle might be alive to
    // another domain
    while (particles.alive(pid) && num_steps < d_max_steps)
    {
        // calculate distance to collision in mean-free-paths
        dist_mfp = -log(particles.ran(pid));

        // while we are hitting boundaries, continue to transport until we get
        // to the collision site
        particles.set_event(pid,profugus::events::BOUNDARY);

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
        while (particles.event(pid) == profugus::events::BOUNDARY)
        {
            DEVICE_CHECK(particles.alive(pid));
            DEVICE_CHECK(dist_mfp > 0.0);

            // total interaction cross section
            xs_tot = d_physics->total(profugus::physics::TOTAL, pid, particles);
            DEVICE_CHECK(xs_tot >= 0.0);

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
            DEVICE_CHECK(step.tag() < profugus::events::END_EVENT);
            particles.set_event(pid,
                static_cast<profugus::events::Event>(step.tag()));

            // path-length tallies (the actual movement of the particle will
            // take place when we process the various events)
            if( d_tallier )
                d_tallier->path_length(step.step(), pid, particles);

            // update the mfp distance travelled
            dist_mfp -= step.step() * xs_tot;

            // process a particle through the geometric boundary
            if (particles.event(pid) == profugus::events::BOUNDARY)
                process_boundary(pid,particles);
        }

        // we are done moving to a problem boundary, now process the
        // collision; the particle is actually moved from its current position
        // in each of these calls

        // process particle at a collision
        if (particles.event(pid) == profugus::events::COLLISION)
            process_collision(pid,particles,step.step());

        // any future events go here ...
        num_steps++;
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

template <class Geometry>
__device__ void
Domain_Transporter<Geometry>::process_boundary(
        int pid, Particle_Vector_t &particles) const
{
    DEVICE_REQUIRE(particles.alive(pid));
    DEVICE_REQUIRE(particles.event(pid) == profugus::events::BOUNDARY);

    // return if not a boundary event

    // move a particle to the surface and process the event through the surface
    d_geometry->move_to_surface(particles.geo_state(pid));

    // get the in/out state of the particle
    int state = d_geometry->boundary_state(particles.geo_state(pid));

    // process the boundary crossing based on the geometric state of the
    // particle
    switch (state)
    {
        case profugus::geometry::OUTSIDE:
        {
            // the particle has left the problem geometry
            particles.set_event(pid,profugus::events::ESCAPE);
            particles.kill(pid);

            DEVICE_ENSURE(particles.event(pid) == profugus::events::ESCAPE);
            break;
        }
        case profugus::geometry::REFLECT:
        {
            // the particle has hit a reflecting surface
            bool reflected = d_geometry->reflect(particles.geo_state(pid));
            DEVICE_CHECK(reflected);

            DEVICE_ENSURE(particles.event(pid) == profugus::events::BOUNDARY);
            break;
        }
        case profugus::geometry::INSIDE:
        {
            // otherwise the particle is at an internal geometry boundary;
            // update the material id of the region the particle has entered
            particles.set_matid(pid,
                                d_geometry->matid(particles.geo_state(pid)));

            // add variance reduction at surface crossings
            if( d_vr )
                d_vr->post_surface(pid, particles);

            DEVICE_ENSURE(particles.event(pid) == profugus::events::BOUNDARY);
            break;

        }
        default:
            DEVICE_CHECK(0);
    }
}

//---------------------------------------------------------------------------//

template <class Geometry>
__device__ void
Domain_Transporter<Geometry>::process_collision(int                pid,
                                                Particle_Vector_t &particles,
                                                double             dist) const
{
    DEVICE_REQUIRE(particles.event(pid) == profugus::events::COLLISION);

    // move the particle to the collision site
    d_geometry->move_to_point(dist, particles.geo_state(pid));

    // sample fission sites
    if (d_sample_fission_sites)
    {
        DEVICE_CHECK(d_fission_sites);
        DEVICE_CHECK(d_keff > 0.0);
        int num_sites = d_physics->sample_fission_site(pid, particles, d_keff);

        if( num_sites > 0 )
        {
            // Create fission site
            Fission_Site site;
            site.m = particles.matid(pid);
            site.r = d_geometry->position(particles.geo_state(pid));

            // Get index of first site to be written using an atomic
            int offset = atomicAdd(d_num_fission_sites,num_sites);
            
            // Don't write past end of allocated array
            num_sites = min(num_sites,d_max_fission_sites-offset);

            // Add fission sites
            for( int i = 0; i < num_sites; ++i )
                d_fission_sites[offset+i] = site;
        }
    }

    // use the physics package to process the collision
    d_physics->collide(pid,particles);

    // apply weight windows
    if( d_vr )
        d_vr->post_collision(pid,particles);
}

} // end namespace cuda_mc

#endif // cuda_mc_Domain_Transporter_i_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.i.cuh
//---------------------------------------------------------------------------//
