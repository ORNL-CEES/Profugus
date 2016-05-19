//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Bank.t.cuh
 * \author Stuart Slattery
 * \date   Friday April 25 16:50:47 2014
 * \brief  Member definitions of class Bank.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Bank_t_cuh
#define cuda_mc_Bank_t_cuh

#include "Bank.hh"

#include "cuda_utils/Memory.cuh"
#include "cuda_utils/Hardware.hh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
template<class Geometry>
__global__
void pop_kernel( const Particle<Geometry>* bank,
		 const std::size_t num_particle,
		 Particle_Vector<Geometry>* particles )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = particles->event_lower_bound( events::DEAD );

    if ( idx < num_particle )
    {

	// Get the particle index.
	int pidx = idx + start_idx;

	// Set the particle data.
	particles->set_wt( pidx, bank[idx].wt() );
	particles->set_group( pidx, bank[idx].group() );
	particles->set_matid( pidx, bank[idx].matid() );
	particles->set_event( pidx, bank[idx].event() );
	particles->geo_state( pidx ) =  bank[idx].geo_state();
	particles->set_batch( pidx, bank[idx].batch() );
        particles->set_dist_mfp( pidx, -std::log(particles->ran(pidx)) );
	particles->live( pidx );
    }    
}

//---------------------------------------------------------------------------//
// HOST API
//---------------------------------------------------------------------------//
// Emit the topmost particles from the stack into empty spots in a vector
template <class Geometry>
void Bank<Geometry>::pop( 
    cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles )
{
    REQUIRE(!empty());
    REQUIRE(!d_particles.empty());
    REQUIRE(!d_count.empty());

    // Get the empty block in the vector.
    int num_particle = particles.get_host_ptr()->get_event_size( events::DEAD );

    // Get the number of particles to pop from the bank.
    int num_to_pop = std::min( num_particle, d_total );

    // Unroll the delayed stack and build a host vector of particles to pass
    // to the device.
    std::vector<Particle_t> host_bank( num_to_pop );
    for ( int i = 0; i < num_to_pop; ++i )
    {
	// Copy the back of the statck into the vector.
	host_bank[i].set_wt( d_particles.back()->wt() );
	host_bank[i].set_group( d_particles.back()->group() );
	host_bank[i].set_matid( d_particles.back()->matid() );
	host_bank[i].set_batch( d_particles.back()->batch() );
	host_bank[i].geo_state() = d_particles.back()->geo_state();
	host_bank[i].set_event( events::TAKE_STEP );

	// Pop the back of the stack either by decrementing the count or
	// actually popping if zero.
	if (d_count.back() > 1)
	{
	    --d_count.back();
	}
	else
	{
	    d_particles.pop_back();
	    d_count.pop_back();
	}
    }

    // Update the running total
    d_total -= num_to_pop;

    // Move the host particles to the device.
    Particle_t* device_bank;
    cuda::memory::Malloc( device_bank, num_to_pop );
    cuda::memory::Copy_To_Device( device_bank, host_bank.data(), num_to_pop );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = num_to_pop / threads_per_block;
    if ( num_to_pop % threads_per_block > 0 ) ++num_blocks;

    // Copy the particles into the vector
    pop_kernel<<<num_blocks,threads_per_block>>>( device_bank,
						  num_to_pop,
						  particles.get_device_ptr() );

    // Free the device bank.
    cuda::memory::Free( device_bank );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Bank_t_cuh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Bank.t.cuh
//---------------------------------------------------------------------------//
