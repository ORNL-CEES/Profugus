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
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Launch_Args.t.cuh"
#include "harness/Diagnostics.hh"
#include "utils/String_Functions.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "Source_Transporter.hh"

#include "Particle_Vector.cuh"
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
    typedef Particle_Vector<Geometry>       Particle_Vector_t;

    // Constructor
    Transport_Functor( const Transporter_t     *trans,
                             Particle_Vector_t  particles,
                       const int               *inds )
        : d_transporter(trans)
        , d_particles(particles)
        , d_inds(inds)
    {
    }

    // Operator apply for functor (transport 1 particle)
    __device__ void operator()( std::size_t tid )
    {
        // Get particle index
        int pid = d_inds[tid];
        DEVICE_CHECK(d_particles.alive(pid));

        // transport the particle through this (replicated) domain
        d_transporter->transport(pid,d_particles);
    }

  private:

    // On-device pointers
    const Transporter_t *d_transporter;
    Particle_Vector_t    d_particles;
    const int           *d_inds;
};

template <class Geometry>
struct GetAlive
{
    typedef Particle_Vector<Geometry> Particle_Vector_t;

    GetAlive(const Particle_Vector_t  particles,
             const int               *indirection,
                   int               *alive)
        : d_particles(particles)
        , d_indirection(indirection)
        , d_alive(alive)
    {
    }

    __device__ void operator()(std::size_t tid)
    {
        int pid = d_indirection[tid];
        if (d_particles.alive(pid))
            d_alive[tid] = 1;
        else
            d_alive[tid] = INT_MAX;
    }

  private:

    const Particle_Vector_t  d_particles;
    const int               *d_indirection;
    int                     *d_alive;

};

template <class Geometry>
struct GetMatid
{
    typedef Particle_Vector<Geometry> Particle_Vector_t;

    GetMatid(const Particle_Vector_t particles,
             const int              *indirection,
                   int              *matids)
        : d_particles(particles)
        , d_indirection(indirection)
        , d_matids(matids)
    {
    }

    __device__ void operator()(std::size_t tid)
    {
        int pid = d_indirection[tid];
        if (d_particles.alive(pid))
            d_matids[tid] = d_particles.matid(pid);
        else
            d_matids[tid] = INT_MAX;
    }

  private:

    const Particle_Vector_t  d_particles;
    const int      *d_indirection;
    int            *d_matids;

};

template <class Geometry>
struct GetGroup
{
    typedef Particle_Vector<Geometry> Particle_Vector_t;

    GetGroup(const Particle_Vector_t particles,
             const int              *indirection,
                   int              *groups)
        : d_particles(particles)
        , d_indirection(indirection)
        , d_groups(groups)
    {
    }

    __device__ void operator()(std::size_t tid)
    {
        int pid = d_indirection[tid];
        if (d_particles.alive(pid))
            d_groups[tid] = d_particles.group(pid);
        else
            d_groups[tid] = INT_MAX;
    }

  private:

    const Particle_Vector_t  d_particles;
    const int               *d_indirection;
    int                     *d_groups;

};

template <class Geometry>
struct GetCell
{
    typedef Particle_Vector<Geometry> Particle_Vector_t;

    GetCell(const Particle_Vector_t   particles,
            const Geometry           *geometry,
            const int                *indirection,
                  int                *cells)
        : d_particles(particles)
        , d_geometry(geometry)
        , d_indirection(indirection)
        , d_cells(cells)
    {
    }

    __device__ void operator()(std::size_t tid)
    {
        int pid = d_indirection[tid];
        if (d_particles.alive(pid))
            d_cells[tid] = d_geometry->cell(d_particles.geo_states(),pid);
        else
            d_cells[tid] = INT_MAX;
    }

  private:

    const Particle_Vector_t  d_particles;
    const Geometry           *d_geometry;
    const int                *d_indirection;
    int                      *d_cells;
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

    d_block_size = db->get("block_size",256);

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

    std::string verb = profugus::to_lower(db->get<std::string>("verbosity",
        std::string("none")));
    if (verb == "none")
        d_verbosity = NONE;
    else if (verb == "low")
        d_verbosity = LOW;
    else if (verb == "medium")
        d_verbosity = MEDIUM;
    else if (verb == "high")
        d_verbosity = HIGH;
    else
        INSIST(false,"Invalid verbosity.");

    d_rng_control = std::make_shared<RNG_Control>(seed);

    // Build domain transporter
    d_transporter = std::make_shared<Transporter_DMM_t>(
        db, d_geometry, d_physics, d_vr );

    // Build particle vector
    d_particle_vec = std::make_shared<Particle_Vector_DMM_t>();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set new tallier
 */
template <class Geometry>
void Source_Transporter<Geometry>::set(SP_Tallier_DMM tallier)
{
    REQUIRE(tallier);
    d_tallier = tallier;
    
    // Don't build on-device tallier here, wait until just before kernel launch
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
    // No global barriers in here.
    // Introducing sync points here will be detrimental to efficiency
    //  of particle batching
    REQUIRE(source);

    SCOPED_TIMER("MC::Source_Transporter.solve");

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> diff;

    // Build Tallier and Domain_Transporter
    auto sdp_tallier = cuda::shared_device_ptr<Tallier<Geometry>>(
        d_tallier->device_instance());
    d_transporter->set(sdp_tallier);
    auto sdp_transporter =
        cuda::shared_device_ptr<Domain_Transporter<Geometry>>(
            d_transporter->device_instance());

    // Get number of particles in source
    size_type num_particles      = source->num_batch();
    size_type num_particles_left = num_particles;

    if (d_verbosity >= LOW)
    {
        std::cout << "Starting solve with " << num_particles
            << " particles on node " << profugus::node() << std::endl;
    }

    thrust::device_vector<int> event(num_particles);
    thrust::device_vector<int> indirection(num_particles);
    thrust::counting_iterator<int> cnt(0);
    thrust::copy(cnt,cnt+num_particles,indirection.begin());

    // Get source particles
    start = std::chrono::high_resolution_clock::now();
    Source_Provider<Geometry> provider;
    provider.get_particles(source, d_rng_control, d_particle_vec, indirection);
    cudaDeviceSynchronize();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    d_source_time += diff.count();

    // Loop until all particles have died
    while (num_particles_left>0)
    {
        // Build launch args
        cuda::Launch_Args<cuda::arch::Device> launch_args;
        launch_args.set_block_size(d_block_size);
        launch_args.set_num_elements(num_particles_left);

        // Build and execute kernel
        cudaDeviceSynchronize();
        start = std::chrono::high_resolution_clock::now();
        Transport_Functor<Geometry> f(
                sdp_transporter.get_device_ptr(), 
                d_particle_vec->device_instance(),
                indirection.data().get() );
        cuda::parallel_launch( f, launch_args );
        REQUIRE(cudaSuccess == cudaGetLastError());
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        d_transport_time += diff.count();

        // Sort particles to put active particles up front
        start = std::chrono::high_resolution_clock::now();
        REQUIRE(cudaSuccess == cudaGetLastError());

        switch (d_sort_type)
        {
            case ALIVE:
            {
                GetAlive<Geometry> alive_func(
                        d_particle_vec->device_instance(),
                        indirection.data().get(),
                        event.data().get());
                cuda::parallel_launch( alive_func, launch_args );
                break;
            }
            case MATID:
            {
                GetMatid<Geometry> matid_func(
                        d_particle_vec->device_instance(),
                        indirection.data().get(),
                        event.data().get());
                cuda::parallel_launch( matid_func, launch_args );
                break;
            }
            case GROUP:
            {
                GetGroup<Geometry> group_func(
                        d_particle_vec->device_instance(),
                        indirection.data().get(),
                        event.data().get());
                cuda::parallel_launch( group_func, launch_args );
                break;
            }
            case CELL:
            {
                GetCell<Geometry> cell_func(
                        d_particle_vec->device_instance(),
                        d_geometry.get_device_ptr(),
                        indirection.data().get(),
                        event.data().get());
                cuda::parallel_launch( cell_func, launch_args );
                break;
            }
        }
        cudaDeviceSynchronize();

        thrust::stable_sort_by_key(event.begin(),
                                   event.begin()+num_particles_left,
                                   indirection.begin());
        auto itr = thrust::find( event.begin(),
                event.begin() + num_particles_left, INT_MAX);
        num_particles_left = itr - event.begin();

        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        d_sort_time += diff.count();

        if (d_verbosity >= HIGH)
        {
            std::cout << num_particles_left << " still alive on node "
                << profugus::node() << std::endl;
        }
    }

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
    d_transporter->set(fis_sites, keff);
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
    return d_transporter->num_sampled_fission_sites();
}


} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.t.cuh
//---------------------------------------------------------------------------//
