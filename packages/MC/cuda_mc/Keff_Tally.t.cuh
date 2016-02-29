//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Keff_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Wed May 14 13:29:40 2014
 * \brief  Keff_Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Keff_Tally_t_hh
#define mc_Keff_Tally_t_hh

#include "Keff_Tally.hh"

#include "comm/global.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Memory.cuh"
#include "cuda_utils/Hardware.hh"
#include "Definitions.hh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Tally accumulate kernel.
template<class Geometry>
__global__ void accumulate_kernel( const Physics<Geometry>* phyiscs,
				   const Particle_Vector<Geometry>* particles,
				   const std::size_t collision_start,
				   const std::size_t num_collision,
				   const std::size_t boundary_start,
				   const std::size_t num_boundary,
				   double* keff )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( (idx >= collision_start && idx < collision_start + num_collision) ||
	 (idx >= boundary_start && idx < boundary_start + num_boundary) )
    {
	keff[idx] = particles->weight(idx) * particles->step(idx) *
		    phyiscs->total( physics::NU_FISSION,
				    particles->matid(idx),
				    particles->group(idx) );
    }
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Kcode solver should construct this with initial keff estimate.
 */
template <class Geometry>
Keff_Tally<Geometry>::Keff_Tally(
    const double keff_init,
    const cuda::Shared_Device_Ptr<Physics_t>& physics )
    : d_physics( physics )
    , d_keff_cycle(keff_init)
{
    REQUIRE(physics);

    // reset tally
    reset();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Geometry>
Keff_Tally<Geometry>::~Keff_Tally()
{
    if ( nullptr != d_keff_device )
    {
	cuda::memory::Free( d_keff_device );
    }
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate average keff over active cycles.
 *
 * This should really only be called after begin_active_cycles; but either way
 * it will return the average keff since either initialization of this tally or
 * the last call to begin_active_cycles.
 *
 * This is only meaningful after one complete cycle. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::mean() const
{
    if (d_cycle < 1)
        return -1.;

    double keff = d_keff_sum / d_cycle;

    ENSURE(keff >= 0.);
    return keff;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate variance of average keff estimate over active cycles.
 *
 * This is only meaningful after two complete cycles. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::variance() const
{
    if (d_cycle < 2)
        return d_keff_sum * d_keff_sum;

    double var = (d_keff_sum_sq - d_keff_sum * d_keff_sum / d_cycle)
                    / (d_cycle * (d_cycle - 1));

    ENSURE(var >= 0.);
    return var;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Track particle and do tallying.
 *
 * This only uses the particle's weight. The accumulated tally is
 * \f[
   k_{\mbox{eff}} = wl\nu\sigma_{\mbox{f}}
 * \f]
 * where \f$l\f$ is the step-length.
 */
template <class Geometry>
void Keff_Tally<Geometry>::accumulate(
    cuda::Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    // Lazy allocate the keff work vector.
    int vector_size = particle.get_host_ptr()->size();
    if ( nullptr == d_keff_device )
    {
	d_keff_host.resize( vector_size );
	cuda::memory::Malloc( d_keff_device, vector_size );
    }
    
    // Get the particles that just had a collision.
    std::size_t collision_start = 0;
    std::size_t num_collision = 0;
    particles.get_host_ptr()->get_event_particles( events::COLLISION,
						   collision_start,
						   num_collision );

    // Get the particles that just hit a boundary.
    std::size_t boundary_start = 0;
    std::size_t num_boundary = 0;
    particles.get_host_ptr()->get_event_particles( events::BOUNDARY,
						   boundary_start,
						   num_boundary );

    // Calculate the launch parameters.
    std::size_t num_particle = num_collision + num_boundary;
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = num_particle / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Process the tallies.
    accumulate_kernel<<<num_blocks,threads_per_block>>>(
	d_physics.get_device_ptr(),
	particles.get_device_ptr(),
	collision_start,
	num_collision,
	boundary_start,
	num_boundary,
	d_keff_device );

    // Pull tallies off the device and add them to cycle tally.
    cuda::memory::Copy_To_Host( 
	d_keff_host.data(), d_keff_device, vector_size );
    for ( int n = collision_start; n < collision_start + num_collision; ++n )
    {
	d_keff_cycle += d_keff_host[n];
    }
    for ( int n = boundary_start; n < boundary_start + num_boundary; ++n )
    {
	d_keff_cycle += d_keff_host[n];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin active cycles in a kcode calculation.
 *
 * This resets the accumulated keff statistics.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_active_cycles()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin a new keff cycle
 *
 * This clears the current accumulated path lengths.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_cycle()
{
    d_keff_cycle = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Mark the end of a keff cycle
 *
 * This performs a global sum across processors, because the provided
 * num_particles is the total number over all domains (blocks plus sets).
 *
 * We accumulate the sum and sum-of-squares of the keff so that we can calculate
 * averages and variances.
 */
template <class Geometry>
void Keff_Tally<Geometry>::end_cycle(double num_particles)
{
    REQUIRE(num_particles > 0.);

    // Keff estimate is total nu-sigma-f reaction rate / num particles
    d_keff_cycle /= num_particles;

    // Do a global sum (since num_particles is global)
    profugus::global_sum(d_keff_cycle);

    // Accumulate first and second moments of cycles, as well as counter
    ++d_cycle;
    d_keff_sum    += d_keff_cycle;
    d_keff_sum_sq += d_keff_cycle * d_keff_cycle;

    // Store keff estimate
    d_all_keff.push_back(d_keff_cycle);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clear accumulated keff values
 *
 * This will still preserve the latest keff value so that existing fission sites
 * will be sampled correctly.
 */
template <class Geometry>
void Keff_Tally<Geometry>::reset()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;

    d_all_keff.clear();
}

} // end namespace cuda_profugus

#endif // mc_Keff_Tally_t_hh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.t.hh
//---------------------------------------------------------------------------//
