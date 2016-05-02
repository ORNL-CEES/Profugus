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
__device__ auto
Uniform_Source<Geometry>::get_particle(std::size_t  tid,
                                       RNG_State_t *rng) const -> Particle_t 
{
    REQUIRE(d_geo_shape);

    // make a particle
    Particle_t p;

    // Set 
    p.set_rng(rng);

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
    b_geometry->initialize(r, omega, p.geo_state());

    // get the material id
    matid = b_geometry->matid(p.geo_state());

    // initialize the physics state by manually sampling the group
    int group = sampler::sample_discrete_CDF(
        d_erg_cdf.size(), d_erg_cdf.data(), p.ran());
    CHECK(group < d_num_groups);
    p.set_group(group);

    // set the material id in the particle
    p.set_matid(matid);

    // set particle weight
    p.set_wt(wt());

    // make particle alive
    p.live();

    ENSURE(p.matid() == matid);
    return p;
}

} // end namespace cuda_mc

#endif // cuda_mc_Uniform_Source_i_cuh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.i.cuh
//---------------------------------------------------------------------------//
