#include <vector>

#include "../Geometry.hh"
#include "../Particle.hh"

void loop_over_particles(acc::Particle *particles);
void ray_trace(acc::Geometry &geometry, int num_rays, int num_steps,
               std::vector<double> &tallies);
