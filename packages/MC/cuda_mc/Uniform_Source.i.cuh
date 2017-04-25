//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.i.cuh
 * \author Steven Hamilton
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source inline member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_i_cuh
#define cuda_mc_Uniform_Source_i_cuh

#include <numeric>

#include "Teuchos_Array.hpp"

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "Sampler.cuh"
#include "Uniform_Source.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Get a particle from the source.
 */
template <class Geometry>
__device__ void Uniform_Source<Geometry>::build_particle(
        int                pid,
        RNG_State_t       *rng,
        Particle_Vector_t &particles) const
{
    DEVICE_REQUIRE(d_geo_shape);

    // Set rng
    particles.set_rng(pid,rng);

    // material id
    int matid = 0;

    // particle position and direction
    cuda_utils::Space_Vector r, omega;

    // sample the angle isotropically
    sampler::sample_isotropic(omega, rng);

    // sample the geometry shape-->we should not get here if there are no
    // particles on this domain
    r = d_geo_shape->sample(rng);

    // intialize the geometry state
    d_geometry->initialize(r, omega, particles.geo_states(),pid);

    // get the material id
    matid = d_geometry->matid(particles.geo_states(),pid);

    // initialize the physics state by manually sampling the group
    int group = sampler::sample_discrete_CDF(
        d_erg_cdf.size(), d_erg_cdf.begin(), particles.ran(pid));
    DEVICE_CHECK(group < d_erg_cdf.size());
    particles.set_group(pid,group);

    // set the material id in the particle
    particles.set_matid(pid,matid);

    // set particle weight
    particles.set_wt(pid,wt());

    // make particle alive
    particles.live(pid);

    DEVICE_ENSURE(particles.matid(pid) == matid);
}

} // end namespace cuda_mc

#endif // cuda_mc_Uniform_Source_i_cuh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.i.cuh
//---------------------------------------------------------------------------//
