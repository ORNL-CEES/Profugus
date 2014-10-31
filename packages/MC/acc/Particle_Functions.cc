//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Particle_Functions.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 28 12:41:23 2014
 * \brief  Particle service functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Particle.hh"

namespace acc
{

//! Hardwired block size.
const def::size_type gang_size   = 60;
const def::size_type vector_size = 128;
const def::size_type rn_per_step = 10;

// Random numbers.
std::vector<double> rnd_numbers;

//---------------------------------------------------------------------------//
/*!
 * \brief Set the number of particles to run concurrently.
 */
void set_size(Vec_Particles &particles,
              int            num_steps)
{
    particles.resize(gang_size * vector_size);
    rnd_numbers.resize(num_steps * rn_per_step);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Load a vector of source particles.
 */
void load_source(profugus::Source &source,
                 Vec_Particles    &particles)
{
    // determine the number of particles in this work block
    def::size_type num_particles = std::min(source.num_to_transport(),
                                            gang_size * vector_size);

    // load the particles from the source
    for (def::size_type p = 0; p < num_particles; ++p)
    {
        // get the particle from the source
        auto particle = source.get_particle();

        // add data to the flat particle
        particles[p].matid = particle->matid();
        particles[p].group = particle->group();
        particles[p].wt    = particle->wt();

        // !!!!!!finish state initialization later!!!!!
    }

    // get a random number generator from the source
    auto rng = source.rng_control().rng();

    // fill up the random numbers
    for (auto &r : rnd_numbers)
    {
        r = rng.ran();
    }
}

} // end namespace acc

//---------------------------------------------------------------------------//
//                 end of Particle_Functions.cc
//---------------------------------------------------------------------------//
