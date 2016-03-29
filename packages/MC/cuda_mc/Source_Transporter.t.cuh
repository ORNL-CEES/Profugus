//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Transporter.t.cuh
 * \author Thomas M. Evans
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_Transporter_t_cuh
#define cuda_mc_Source_Transporter_t_cuh

#include <iomanip>
#include <iostream>
#include <cmath>
#include <thrust/device_vector.h>

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Launch_Args.t.cuh"
#include "harness/Diagnostics.hh"
#include "utils/String_Functions.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "Source_Transporter.hh"

#include "Particle.cuh"
#include "Physics.cuh"
#include "Domain_Transporter.cuh"
#include "VR_Roulette.cuh"
#include "Tallier.cuh"

namespace cuda_mc
{

// Functor to transport source particles
template <class Geometry>
class Transport_Functor
{
  public:

    typedef Domain_Transporter<Geometry>    Transporter_t;
    typedef Tallier<Geometry>               Tallier_t;
    typedef Particle<Geometry>              Particle_t;

    // Constructor
    Transport_Functor( const Transporter_t *trans,
                             Particle_t    *particles )
        : d_transporter( trans )
        , d_particles( particles )
    {
    }

    // Operator apply for functor (transport 1 particle)
    __device__ void operator()( std::size_t tid ) const
    {
        // Get particle from source
        auto &p = d_particles[tid];
        CHECK( p.alive() );

        // transport the particle through this (replicated) domain
        d_transporter->transport(p);
        CHECK(!p.alive());
    }

  private:

    // On-device pointers
    const Transporter_t *d_transporter;
    Particle_t          *d_particles;
};

// Functor to initialize RNG
class RNG_Init
{
  public:

    RNG_Init( RNG_State_t *rngs)
      : d_rngs(rngs)
    {}

    __device__ void operator()( std::size_t tid ) const
    {
        curand_init(1234,tid,0,d_rngs+tid);
    }

  private:

    RNG_State_t *d_rngs;
};
        

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor/
 */
template <class Geometry>
Source_Transporter<Geometry>::Source_Transporter(RCP_Std_DB   db,
                                                 SDP_Geometry geometry,
                                                 SDP_Physics  physics,
                                                 SDP_Tallier  tallier)
    : d_geometry(geometry)
    , d_physics(physics)
    , d_tallier(tallier)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
    REQUIRE(!db.is_null());
    REQUIRE(d_geometry.get_host_ptr());
    REQUIRE(d_geometry.get_device_ptr());
    REQUIRE(d_physics.get_host_ptr());
    REQUIRE(d_physics.get_device_ptr());

    // Build variance reduction
    std::string var = profugus::to_lower(
        db->get<std::string>("variance reduction",std::string("roulette")) );
    if( var == "roulette" )
    {
        d_vr = cuda::shared_device_ptr<VR_Roulette_t>(db);
    }

    // Build domain transporter
    d_transporter = cuda::shared_device_ptr<Transporter_t>(
        d_geometry, d_physics, d_tallier, d_vr );
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 */
template <class Geometry>
template <class Src_Type>
void
Source_Transporter<Geometry>::solve(std::shared_ptr<Src_Type> source) const
{
    REQUIRE(source);

    // barrier at the start
    profugus::global_barrier();

    SCOPED_TIMER("MC::Source_Transporter.solve");

    // Copy class to device
    // Is it better to have this method take an SDP<Src_Type> and push
    // the copy-to-device outside of the Source_Transporter?
    cuda::Shared_Device_Ptr<Src_Type> sdp_source(source);
    ENSURE( sdp_source.get_device_ptr() );

    // Get number of particles in source
    int num_particles = source->num_to_transport();

    // Build launch args
    cuda::Launch_Args<cuda::arch::Device> launch_args;
    launch_args.set_num_elements(num_particles);

    // Initialize RNG
    thrust::device_vector<RNG_State_t> rngs(num_particles);
    RNG_Init init( rngs.data().get() );
    cuda::parallel_launch( init, launch_args );
    cudaDeviceSynchronize();

    // Get source particles
    auto particles = get_particles( sdp_source, rngs );
    CHECK( particles.size() == num_particles );
    cudaDeviceSynchronize();

    // Build and execute kernel
    Transport_Functor<Geometry> f( d_transporter.get_device_ptr(), 
                                   particles.data().get() );
    cuda::parallel_launch( f, launch_args );
    cudaDeviceSynchronize();

    // barrier at the end
    profugus::global_barrier();

    // increment the particle counter
    DIAGNOSTICS_ONE(integers["particles_transported"] += num_particles);

#ifdef REMEMBER_ON
    profugus::global_sum(num_particles);
    ENSURE(num_particles== source->total_num_to_transport());
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set a fission site container and keff for sampling fission sites.
 *
 * Setting a fission site container tells the fixed-source solver to sample
 * fission sites that can be used in an outer k-code calculation.  Thus, this
 * function should be called if the fixed-source solver is used as the inner
 * part of a k-code eigenvalue calculation.  It should be called once per
 * k-code iteration to update the eigenvalue.
 *
 * Fission sites are added to the container, it is \b not emptied.
 */
#if 0
template <class Geometry>
void Source_Transporter<Geometry>::sample_fission_sites(SP_Fission_Sites fis_sites,
                                              double           keff)
{
    // set the transporter with the fission site container and the latest keff
    // iterate
    d_transporter.set(fis_sites, keff);
}
#endif

} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.t.cuh
//---------------------------------------------------------------------------//
