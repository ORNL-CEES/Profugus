#include <cmath>

#include "utils/Constants.hh"
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

void ray_trace(acc::Geometry             &geometry,
               int                        num_rays,
               const std::vector<double> &rnd,
               std::vector<double>       &tallies)
{
    const double *rnd_ptr = &rnd[0];
    int rctr = 0;

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
            r.pos[d] = extent * rnd[rctr++];

            // find ijk position
            r.ijk[d] = r.pos[d] / delta;
        }

        // sample direction
        double costheta = 2.0 * rnd[rctr++] - 1.0;
        double sintheta = std::sqrt(1.0 - costheta * costheta);
        double phi      = profugus::constants::two_pi * rnd[rctr++];

        r.dir[0] = sintheta * std::cos(phi);
        r.dir[1] = sintheta * std::sin(phi);
        r.dir[2] = costheta;
    }

    // get pointer to rays
    acc::Geometry_State *ray_ptr = &rays[0];

    // get pointer to tallies
    double *tally = &tallies[0];
    int nc        = geometry.num_cells();

#pragma acc parallel loop present(geometry) copyin(ray_ptr[0:num_rays]) \
    copy(tally[0:nc])
    {
        for (int r = 0; r < num_rays; ++r)
        {
            // get reference reference to ray
            acc::Geometry_State &ray = ray_ptr[r];

            // loop over steps for each ray
#pragma acc loop seq
            for (int s = 0; s < 1000; ++s)
            {
                double dbnd = geometry.distance_to_boundary(ray);

                // get the cell index
                int cell = geometry.cell(ray);

                // tally the pathlength
                tally[cell] += dbnd;

                // move the ray to the next surface
                geometry.move_to_surface(ray);

                // reflect the particle
                if (geometry.boundary_state(ray) == profugus::geometry::REFLECT)
                {
                    geometry.reflect(ray);
                }
            }
        }
    }
}
