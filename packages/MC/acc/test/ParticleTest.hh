#include <vector>

#include "../Geometry.hh"
#include "../Particle.hh"

void loop_over_particles(acc::Particle *particles);
void ray_trace(acc::Geometry &geometry, int num_rays,
               const std::vector<double> &rnd, std::vector<double> &tallies);
