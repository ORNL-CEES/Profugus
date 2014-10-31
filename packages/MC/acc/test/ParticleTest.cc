#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "utils/Constants.hh"
#include "../RNG.hh"
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

//---------------------------------------------------------------------------//

void ray_trace(acc::Geometry       &geometry,
               int                  num_rays,
               int                  num_steps,
               std::vector<double> &tallies)
{
#ifndef _OPENACC
    using std::log;
#endif

    acc::RNG rng(32423);

    // seeds
    int ctr = 0;
    std::vector<long> seeds(num_rays, 200);

    // make a vector of states
    std::vector<acc::Geometry_State> rays(num_rays);
    for (auto &r : rays)
    {
        // sample a position in the geometry
        for (int d = 0; d < 3; ++d)
        {
            int N         = geometry.num_cells(d);
            double extent = geometry.extents(d)[N];
            double delta  = extent / N;

            // sample position
            r.pos[d] = extent * rng.ran();

            // find ijk position
            r.ijk[d] = r.pos[d] / delta;
        }

        // sample direction
        double costheta = 2.0 * rng.ran() - 1.0;
        double sintheta = std::sqrt(1.0 - costheta * costheta);
        double phi      = profugus::constants::two_pi * rng.ran();

        r.dir[0] = sintheta * std::cos(phi);
        r.dir[1] = sintheta * std::sin(phi);
        r.dir[2] = costheta;

        seeds[ctr] += ctr;
        ++ctr;
    }

    // get pointer to rays
    acc::Geometry_State *ray_ptr = &rays[0];
    long *seeds_ptr              = &seeds[0];

    // get pointer to tallies
    double *tally  = &tallies[0];
    int nc         = geometry.num_cells();
    double keff    = 0.0;

    // weights
    std::vector<double> wts(num_rays, 1.0);
    double *wts_ptr = &wts[0];

#pragma acc parallel loop present(geometry) copyin(ray_ptr[0:num_rays]) \
    copy(tally[0:nc]) present(rng) pcopyin(seeds_ptr[0:num_rays])       \
    pcopyin(wts_ptr[0:num_rays])
    {
        for (int r = 0; r < num_rays; ++r)
        {
            // get reference reference to ray
            acc::Geometry_State &ray = ray_ptr[r];

            // k-tally
            double k = 0.0;

            // step-length
            double step = 0.0;
            int type    = 0;

            // loop over steps for each ray
#pragma acc loop seq
            for (int s = 0; s < num_steps; ++s)
            {
                // sample distance to collision
                double dcol = (-5.0 * log(rng.ran(seeds_ptr[r])));

                // distance to boundary
                double dbnd = geometry.distance_to_boundary(ray);

                // determine step
                step = dbnd;
                type = 0;
                if (dcol < dbnd)
                {
                     step = dcol;
                     type = 1;
                }

                // update distance to collision
                //dcol -= step;

                // global tally
                k += step * wts_ptr[r];

                if (type == 1)
                {
                    wts_ptr[r] = 0.0;
                }

                // move the ray to the next surface
                geometry.move_to_surface(ray);

                // reflect the particle
                if (geometry.boundary_state(ray) == profugus::geometry::REFLECT)
                {
                    geometry.reflect(ray);
                }
            }

            keff += k;
        }
    }

    tallies[0] = keff;
}
