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
__device__ auto Uniform_Source<Geometry>::get_particle(RNG_t &rng) -> Particle_t 
{
    REQUIRE(d_wt > 0.0);
    REQUIRE(d_geo_shape);

    // make a particle
    Particle_t p;

    // Set 
    p.set_rng(rng);

    // material id
    int matid = 0;

    // particle position and direction
    cuda::Space_Vector r, omega;

    // sample the angle isotropically
    Base::sample_angle(omega, rng);

    // sample the geometry shape-->we should not get here if there are no
    // particles on this domain
    r = d_geo_shape->sample(rng);

    // intialize the geometry state
    b_geometry->initialize(r, omega, p.geo_state());

    // get the material id
    matid = b_geometry->matid(p.geo_state());

    if( blockIdx.x==0 && threadIdx.x==0 )
    {
        printf("Energy cdf: ");
        double *dat = d_erg_cdf.data();
        for( int i = 0; i < d_erg_cdf.size(); ++i )
            printf("%e, ",dat[i]);
        printf("\n");
    }

    // initialize the physics state by manually sampling the group
    int group = sampler::sample_discrete_CDF(
        d_erg_cdf.size(), d_erg_cdf.data(), p.ran());
    CHECK(group < d_num_groups);
    p.set_group(group);

    // set the material id in the particle
    p.set_matid(matid);

    // set particle weight
    p.set_wt(d_wt);

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
