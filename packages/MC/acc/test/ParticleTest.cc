#include "ParticleTest.hh"

void loop_over_particles(acc::Particle *particles)
{
#pragma acc kernels copyout(particles[0:48])
    for (int n = 0; n < 48; ++n)
    {
        // get a particle
        acc::Particle &p = particles[n];
    }
}
