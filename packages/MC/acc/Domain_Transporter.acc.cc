//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Domain_Transporter.acc.cc
 * \author Thomas M. Evans
 * \date   Sat Nov 01 10:56:35 2014
 * \brief  Domain_Transporter OpenACC member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Domain_Transporter.hh"

namespace acc
{

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Domain_Transporter::~Domain_Transporter()
{
}

//---------------------------------------------------------------------------//
// ACCELERATED FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Transport a vector of ACC-device particles.
 *
 * This function is designed to be called from the CPU, not the device!
 */
void Domain_Transporter::transport(Vec_ACC_Particles &particles)
{
    // number of particles to loop over
    int num_particles = particles.size();

    // get pointers/references to objects
    Geometry &geometry      = *d_geometry;
    Particle *particles_ptr = &particles[0];

#pragma acc parallel loop present(geometry) present(d_rng)                      \
    copyin(particles_ptr[0:num_particles]
    {
        for (int n = 0; n < num_particles; ++n)
        {
            // get a reference to the particle
            Particle &p = particles_ptr[n];
        }
    }
}

} // end namespace acc

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.acc.cc
//---------------------------------------------------------------------------//
