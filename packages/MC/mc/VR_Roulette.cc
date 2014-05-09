//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/VR_Roulette.cc
 * \author Thomas M. Evans
 * \date   Fri May 09 13:09:37 2014
 * \brief  VR_Roulette member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "VR_Roulette.hh"
#include "harness/Diagnostics.hh"
#include "Definitions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
VR_Roulette::VR_Roulette(RCP_Std_DB db)
    : Base()
{
    Require (!db.is_null());

    // Get default weight cutoffs and survival
    d_Wc = db->get("weight_cutoff", 0.25);
    d_Ws = db->get("weight_survival", 2.0 * d_Wc);

    b_splitting = false;

    Validate(d_Wc >= 0, "Weight cutoff must be nonnegative "
            "(got weight_cutoff=" << d_Wc << ")");

    Validate(d_Ws > 0 ? d_Ws > d_Wc : d_Wc == 0.,
            "Survival weight must be greater than cutoff weight (got "
            "weight_cutoff=" << d_Wc << ", "
            "weight_survival=" << d_Ws << ")");
}

//---------------------------------------------------------------------------//
// VARIANCE REDUCTION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle for weight roulette.
 */
void VR_Roulette::post_collision(Particle_t& particle,
                                 Bank_t&     bank) const
{
    if (!particle.alive())
        return;

    // get the particle weight
    const double orig_weight = particle.wt();

    // if the particle weight is below the cutoff do roulette
    if (orig_weight < d_Wc)
    {
        // the particle should always be alive if it gets here
        Check (particle.alive());
        Check (d_Ws >= d_Wc);

        // calculate survival probablity
        const double survival = orig_weight / d_Ws;
        Check (survival < 1.0);

        // particle survives roulette
        if (particle.rng().ran() < survival)
        {
            // set the new weight of the surviving particle
            particle.set_wt(d_Ws);
            Check (particle.wt() == d_Ws);

            // update the event
            particle.set_event(events::ROULETTE_SURVIVE);

            // add the event to the diagnostics
            DIAGNOSTICS_TWO(integers["roulette_survive"]++);
        }

        // otherwise the particle dies
        else
        {
            // kill the particle
            particle.kill();

            // update the event
            particle.set_event(events::ROULETTE_KILLED);

            // add the event to the diagnostics
            DIAGNOSTICS_TWO(integers["roulette_killed"]++);
        }
    }

    Ensure (particle.wt() >= orig_weight);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.cc
//---------------------------------------------------------------------------//
