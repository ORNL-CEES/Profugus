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
auto Fission_Source<Geometry>::get_particle(
    std::size_t tid, RNG_State_t *rng) const -> Particle_t
{
    DEVICE_REQUIRE(d_wt > 0.0);

    // particle
    Particle_t p;
    DEVICE_CHECK(!p.alive());

    // use the global rng on this domain for the random number generator
    p.set_rng(rng);

    // material id
    int matid = 0;

    // particle position and isotropic direction
    Space_Vector omega;

    // sample the angle isotropically
    sampler::sample_isotropic(omega, rng);

    // if there is a fission site container than get the particle from there;
    // otherwise assume this is an initial source
    if (d_have_sites)
    {
        // get the corresponding element in the site container
        Fission_Site &fs = d_fission_sites[tid];

        // intialize the geometry state
        b_geometry->initialize(fs.r, omega, p.geo_state());

        // get the material id
        matid = b_geometry->matid(p.geo_state());

        // initialize the physics state at the fission site
        bool sampled = d_physics->initialize_fission(fs, p);
        DEVICE_CHECK(sampled);
    }
    else
    {
        Space_Vector r;
        matid = sample_geometry(r, omega, p, rng);
    }

    // set the material id in the particle
    p.set_matid(matid);

    // set particle weight
    p.set_wt(d_wt);

    // make particle alive
    p.live();

    DEVICE_ENSURE(p.matid() == matid);
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample geometry to get a particle.
 */
template <class Geometry>
__device__
int Fission_Source<Geometry>::sample_geometry(Space_Vector       &r,
                                              const Space_Vector &omega,
                                              Particle_t         &p,
                                              RNG_State_t        *rng) const
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
        r[I] = d_width[I] * curand_uniform_double(rng) + d_lower[I];
        r[J] = d_width[J] * curand_uniform_double(rng) + d_lower[J];
        r[K] = d_width[K] * curand_uniform_double(rng) + d_lower[K];

        // intialize the geometry state
        b_geometry->initialize(r, omega, p.geo_state());

        // get the material id
        matid = b_geometry->matid(p.geo_state());

        // try initializing fission here, if it is successful we are
        // finished
        if (d_physics->initialize_fission(matid, p))
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
