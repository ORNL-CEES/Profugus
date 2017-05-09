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
void VR_Roulette<Geometry>::post_collision(int                pid,
                                           Particle_Vector_t* particles) const
{
    if (!particles->alive(pid))
        return;

    // get the particle weight
    const double orig_weight = particles->wt(pid);

    // if the particle weight is below the cutoff do roulette
    if (orig_weight < d_Wc)
    {
        // the particle should always be alive if it gets here
        DEVICE_CHECK(particles->alive(pid));
        DEVICE_CHECK(d_Ws >= d_Wc);

        // calculate survival probablity
        const double survival = orig_weight / d_Ws;
        DEVICE_CHECK(survival < 1.0);

        // particle survives roulette
        if (particles->ran(pid) < survival)
        {
            // set the new weight of the surviving particle
            particles->set_wt(pid,d_Ws);
            DEVICE_CHECK(particles->wt(pid) == d_Ws);

            // update the event
            particles->set_event(pid,profugus::events::ROULETTE_SURVIVE);
        }

        // otherwise the particle dies
        else
        {
            // kill the particle
            particles->kill(pid);

            // update the event
            particles->set_event(pid,profugus::events::ROULETTE_KILLED);
        }
    }

    DEVICE_ENSURE(particles->wt(pid) >= orig_weight);
}

} // end namespace cuda_mc 

#endif // cuda_mc_VR_Roulette_i_cuh

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.i.cuh
//---------------------------------------------------------------------------//
