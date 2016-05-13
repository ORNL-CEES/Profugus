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
#include <chrono>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/count.h>

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
#include "RNG_Control.cuh"
#include "Source_Provider.cuh"

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
    }

  private:

    // On-device pointers
    const Transporter_t *d_transporter;
    Particle_t          *d_particles;
};

template <class Particle>
struct AliveComp
{
    __device__ bool operator()(const Particle &lhs, const Particle &rhs) const
    {
        return lhs.alive() && !rhs.alive();
    }
};

template <class Particle>
struct MatidComp 
{
    __device__ bool operator()(const Particle &lhs, const Particle &rhs) const
    {
        return comp(lhs,rhs) || ((lhs.alive() && rhs.alive()) &&
               (lhs.matid() < rhs.matid()));
    }

    AliveComp<Particle> comp;
};

template <class Particle>
struct GroupComp 
{
    __device__ bool operator()(const Particle &lhs, const Particle &rhs) const
    {
        return comp(lhs,rhs) || ((lhs.alive() && rhs.alive()) &&
               (lhs.group() < rhs.group()));
    }

    AliveComp<Particle> comp;
};

template <class Geometry>
struct CellComp 
{
    typedef Particle<Geometry> Particle_t;
    __device__ bool operator()(const Particle_t &lhs, const Particle_t &rhs) const
    {
        return comp(lhs,rhs) || ((lhs.alive() && rhs.alive()) &&
               (geom->cell(lhs.geo_state()) < geom->cell(rhs.geo_state())));
    }

    AliveComp<Particle_t> comp;
    const Geometry *geom;
};

template <class Particle>
struct IsAlive
{
    __host__ __device__ bool operator()(const Particle &p) const
    {
        return p.alive();
    }
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
    , d_sort_time(0.0)
    , d_transport_time(0.0)
    , d_source_time(0.0)
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

    int seed = db->get("seed",1234);

    std::string sort_type = profugus::to_lower(db->get<std::string>("sort_type",
            std::string("alive")));
    if (sort_type == "alive" )
        d_sort_type = ALIVE;
    else if (sort_type == "matid")
        d_sort_type = MATID;
    else if (sort_type == "group")
        d_sort_type = GROUP;
    else if (sort_type == "cell")
        d_sort_type = CELL;
    else
        INSIST(false,"Invalid sort type.");

    d_rng_control = std::make_shared<RNG_Control>(seed);

    // Build domain transporter
    d_transporter = cuda::shared_device_ptr<Transporter_t>(
        db, d_geometry, d_physics, d_vr );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set new tallier
 */
template <class Geometry>
void Source_Transporter<Geometry>::set(SDP_Tallier tallier)
{
    REQUIRE( tallier.get_host_ptr() );
    REQUIRE( tallier.get_device_ptr() );
    d_tallier = tallier;

    d_transporter.get_host_ptr()->set(d_tallier);
    d_transporter.update_device();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 */
template <class Geometry>
void Source_Transporter<Geometry>::solve(SP_Source source) const
{
    REQUIRE(source);

    // barrier at the start
    profugus::global_barrier();

    SCOPED_TIMER("MC::Source_Transporter.solve");

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> diff;

    // Get source particles
    start = std::chrono::high_resolution_clock::now();
    Source_Provider<Geometry> provider;
    thrust::device_vector<Particle_t> particles;
    provider.get_particles( source, d_rng_control, particles );
    cudaDeviceSynchronize();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    d_source_time += diff.count();

    // Get number of particles in source
    int num_particles = particles.size();

    std::cout << "Starting with " << num_particles << " particles" << std::endl;
    while (num_particles>0)
    {
        // Build launch args
        cuda::Launch_Args<cuda::arch::Device> launch_args;
        launch_args.set_num_elements(num_particles);

        // Build and execute kernel
        cudaDeviceSynchronize();
        start = std::chrono::high_resolution_clock::now();
        Transport_Functor<Geometry> f( d_transporter.get_device_ptr(), 
                                       particles.data().get() );
        cuda::parallel_launch( f, launch_args );
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        d_transport_time += diff.count();

        // Sort particles to put active particles up front
        start = std::chrono::high_resolution_clock::now();

        switch (d_sort_type)
        {
            case ALIVE:
                thrust::sort( particles.begin(),
                              particles.begin() + num_particles,
                              AliveComp<Particle_t>() );
                break;
            case MATID:
                thrust::sort( particles.begin(),
                              particles.begin() + num_particles,
                              MatidComp<Particle_t>() );
                break;
            case GROUP:
                thrust::sort( particles.begin(),
                              particles.begin() + num_particles,
                              GroupComp<Particle_t>() );
                break;
            case CELL:
                CellComp<Geometry> comp;
                comp.geom = d_geometry.get_device_ptr();
                thrust::sort(particles.begin(),
                             particles.begin() + num_particles,
                             comp);
                break;
        }

        // Count particles still alive
        num_particles = thrust::count_if( particles.begin(),
                particles.begin() + num_particles, IsAlive<Particle_t>() );
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        d_sort_time += diff.count();

        std::cout << num_particles << " still alive" << std::endl;

    }

    // barrier at the end
    profugus::global_barrier();

    // increment the particle counter
    DIAGNOSTICS_ONE(integers["particles_transported"] += num_particles);
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
template <class Geometry>
void Source_Transporter<Geometry>::sample_fission_sites(
    SP_Fission_Site_Vec fis_sites, double keff)
{
    // set the transporter with the fission site container and the latest keff
    // iterate
    d_transporter.get_host_ptr()->set(fis_sites, keff);
    d_transporter.update_device();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get number of fission sites sampled during last transport.
 *
 * The size of the fission site vector is not changed during transport
 * (dynamic reallocation inside a kernel is slow and dangerous), so we 
 * need to resize the vector to the number of sites that were actually created.
 * Here we get the number of sampled sites from the Domain_Transporter.
 * Note that this requires updating the host-side copy of the
 * Domain_Transporter because the number of fission sites is accumulated as
 * a member variable during transport.
 *
 */
template <class Geometry>
int Source_Transporter<Geometry>::num_sampled_fission_sites()
{
    d_transporter.update_host();
    int num_fission_sites =
        d_transporter.get_host_ptr()->num_sampled_fission_sites();
    return num_fission_sites;
}


} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.t.cuh
//---------------------------------------------------------------------------//
