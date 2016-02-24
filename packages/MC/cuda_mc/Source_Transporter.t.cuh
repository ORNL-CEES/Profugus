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

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Launch_Args.hh"
#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "Source_Transporter.hh"

#include "Particle.cuh"
#include "Physics.cuh"
#include "Domain_Transporter.cuh"

namespace cuda_mc
{

// Functor to transport source particles
template <class Geometry, class Src_Type>
class Source_Functor
{
  public:

    typedef Domain_Transporter<Geometry>            Transporter_t;
    typedef cuda::Shared_Device_Ptr<Src_Type>       SDP_Source;
    typedef cuda::Shared_Device_Ptr<Transporter_t>  SDP_Transporter;

    Source_Functor( SDP_Source      src,
                    SDP_Transporter trans )
        : d_src( src.get_device_ptr() )
        , d_transporter( trans.get_device_ptr() )
    {
    }

    __device__ void operator()( const std::size_t tid )
    {
        // Get particle from source
        auto p = d_src->get_particle(tid);
        CHECK( p.alive() );

        // Do "source event" tallies on the particle
        //d_tallier->source(p);

        // transport the particle through this (replicated) domain
        d_transporter->transport(p);
        CHECK(!p.alive());

        // indicate completion of particle history
        //d_tallier->end_history();
    }

  private:

    // On-device pointers
    Src_Type        *d_src;
    Transporter_t   *d_transporter;
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
                                                 SDP_Physics  physics)
    : d_geometry(geometry)
    , d_physics(physics)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
    REQUIRE(!db.is_null());
    REQUIRE(d_geometry.get_host_ptr());
    REQUIRE(d_geometry.get_device_ptr());
    REQUIRE(d_physics.get_host_ptr());
    REQUIRE(d_physics.get_device_ptr());

    // set the geometry and physics in the domain transporter
    auto trans_host = std::make_shared<Transporter_t>();
    trans_host->set(d_geometry, d_physics);
    d_transporter = SDP_Transporter( trans_host );
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

    // Build and execute kernel
    Source_Functor<Geometry,Src_Type> f( d_transporter, sdp_source );
    cuda::parallel_launch( f, launch_args );

    // barrier at the end
    profugus::global_barrier();

    // increment the particle counter
    DIAGNOSTICS_ONE(integers["particles_transported"] += num_particles);

#ifdef REMEMBER_ON
    profugus::global_sum(num_particles);
    ENSURE(num_particles== source.total_num_to_transport());
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

//---------------------------------------------------------------------------//
/*!
 * \brief Set the variance reduction.
 */
#if 0
template <class Geometry>
void Source_Transporter<Geometry>::set(SP_Variance_Reduction vr)
{
    REQUIRE(vr);

    // set the variance reduction in the domain transporter and locally
    d_transporter.set(vr);
    d_var_reduction = vr;

    ENSURE(d_var_reduction);
}
#endif

//---------------------------------------------------------------------------//
/*!
 * \brief Set the tally controller
 */
#if 0
template <class Geometry>
void Source_Transporter<Geometry>::set(SP_Tallier tallier)
{
    REQUIRE(tallier);

    // set the tally controller in the domain transporter and locally
    d_transporter.set(tallier);
    d_tallier = tallier;

    ENSURE(d_tallier);
}
#endif

} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.t.cuh
//---------------------------------------------------------------------------//
