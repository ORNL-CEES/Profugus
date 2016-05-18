//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/VR_Roulette.i.cuh
 * \author Thomas M. Evans
 * \date   Fri May 09 13:09:37 2014
 * \brief  VR_Roulette member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_VR_Roulette_i_cuh
#define cuda_mc_VR_Roulette_i_cuh

#include "VR_Roulette.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// VARIANCE REDUCTION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle for weight roulette.
 */
template <class Geometry>
__device__
void VR_Roulette<Geometry>::post_collision(Particle_t& particle) const
{
    if (!particle.alive())
        return;

    // get the particle weight
    const double orig_weight = particle.wt();

    // if the particle weight is below the cutoff do roulette
    if (orig_weight < d_Wc)
    {
        // the particle should always be alive if it gets here
        DEVICE_CHECK(particle.alive());
        DEVICE_CHECK(d_Ws >= d_Wc);

        // calculate survival probablity
        const double survival = orig_weight / d_Ws;
        DEVICE_CHECK(survival < 1.0);

        // particle survives roulette
        if (particle.ran() < survival)
        {
            // set the new weight of the surviving particle
            particle.set_wt(d_Ws);
            DEVICE_CHECK(particle.wt() == d_Ws);

            // update the event
            particle.set_event(profugus::events::ROULETTE_SURVIVE);
        }

        // otherwise the particle dies
        else
        {
            // kill the particle
            particle.kill();

            // update the event
            particle.set_event(profugus::events::ROULETTE_KILLED);
        }
    }

    DEVICE_ENSURE(particle.wt() >= orig_weight);
}

} // end namespace cuda_mc 

#endif // cuda_mc_VR_Roulette_i_cuh

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.i.cuh
//---------------------------------------------------------------------------//
