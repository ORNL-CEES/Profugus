//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.i.cuh
 * \author Steven Hamilton
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source inline member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fission_Source_i_cuh
#define cuda_mc_Fission_Source_i_cuh

#include <algorithm>
#include <numeric>

#include "Teuchos_Array.hpp"

#include "Fission_Source.cuh"
#include "Sampler.cuh"

#include "harness/DBC.hh"
#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "utils/Constants.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Get a particle from the source.
*/
template <class Geometry>
__device__
void Fission_Source<Geometry>::build_particle(
    int pid, Particle_Vector_t *particles) const
{
    DEVICE_REQUIRE(d_wt > 0.0);

    // material id
    int matid = 0;

    // particle position and isotropic direction
    Space_Vector omega;

    // sample the angle isotropically
    sampler::sample_isotropic(omega, particles->rng(pid));

    // if there is a fission site container than get the particle from there;
    // otherwise assume this is an initial source
    if (!is_initial_source())
    {
        // get the corresponding element in the site container
        const Fission_Site &fs = d_fission_sites[cuda::utility::thread_id()];

        // intialize the geometry state
        d_geometry->initialize(fs.r, omega, particles->geo_states(),pid);

        // get the material id
        matid = d_geometry->matid(particles->geo_states(),pid);

        // initialize the physics state at the fission site
        bool sampled = d_physics->initialize_fission(fs, pid, particles);
        DEVICE_CHECK(sampled);
    }
    else
    {
        Space_Vector r;
        matid = sample_geometry(r, omega, pid, particles);
    }

    // set the material id in the particle
    particles->set_matid(pid,matid);

    // set particle weight
    particles->set_wt(pid,d_wt);

    // make particle alive
    particles->live(pid);

    DEVICE_ENSURE(particles->matid(pid) == matid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample geometry to get a particle.
 */
template <class Geometry>
__device__
int Fission_Source<Geometry>::sample_geometry(Space_Vector       &r,
                                              const Space_Vector &omega,
                                              int                 pid,
                                              Particle_Vector_t  *particles) const
{
    using def::I; using def::J; using def::K;

    // sampled complete flag
    bool sampled = false;

    // material id
    int matid = 0;

    // >>> Sample the full geometry

    // sample the geometry until a fission site is found (if there is no
    // fission in a given domain the number of particles on that domain is
    // zero, and we never get here) --> so, fission sampling should always
    // be successful
    while (!sampled)
    {
        // sample a point in the geometry
        r[I] = d_width[I] * particles->ran(pid) + d_lower[I];
        r[J] = d_width[J] * particles->ran(pid) + d_lower[J];
        r[K] = d_width[K] * particles->ran(pid) + d_lower[K];

        // intialize the geometry state
        d_geometry->initialize(r, omega, particles->geo_states(),pid);

        // get the material id
        matid = d_geometry->matid(particles->geo_states(),pid);

        // try initializing fission here, if it is successful we are
        // finished
        if (d_physics->initialize_fission(matid, pid, particles))
        {
            sampled = true;
        }
    }

    return matid;
}

} // end namespace cuda_mc

#endif // cuda_mc_Fission_Source_i_cuh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.i.cuh
//---------------------------------------------------------------------------//
