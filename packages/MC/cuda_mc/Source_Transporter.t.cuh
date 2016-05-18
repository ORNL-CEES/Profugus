//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source_Transporter.t.hh
 * \author Stuart Slattery
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

#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "mc/Global_RNG.hh"
#include "Source_Transporter.hh"

#include <cuda_profiler_api.h>

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor/
 */
template <class Geometry>
Source_Transporter<Geometry>::Source_Transporter(const RCP_Std_DB& db,
                                                 const SDP_Geometry& geometry,
                                                 const SDP_Physics& physics)
    : d_geometry(geometry)
    , d_physics(physics)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
    REQUIRE(!db.is_null());
    REQUIRE(d_geometry);
    REQUIRE(d_physics);

    // set the geometry and physics in the domain transporter
    d_transporter.set(d_geometry, d_physics);

    // set the output frequency for particle transport diagnostics
    d_print_fraction = db->get("cuda_mc_diag_frac", 1.1);

    // Get the particle vector size.
    d_vector_size = db->get("particle_vector_size",10000);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Assign the source.
 */
template <class Geometry>
void Source_Transporter<Geometry>::assign_source(const SP_Source& source)
{
    using std::ceil;

    REQUIRE(source);

    // assign the source
    d_source = source;

    // calculate the frequency of output diagnostics
    d_print_count = ceil(d_source->num_to_transport() * d_print_fraction);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 */
template <class Geometry>
void Source_Transporter<Geometry>::solve()
{
    using std::cout; using std::endl;

    REQUIRE(d_source);

    // barrier at the start
    profugus::global_barrier();

    SCOPED_TIMER("CUDA_MC::Source_Transporter.solve");

    // event counter
    size_type counter = 0;

    // make a particle bank
    cuda::Shared_Device_Ptr<typename Transporter_t::Bank_t> bank;

    // make a particle vector
    auto particles = cuda::shared_device_ptr<Particle_Vector<Geometry> >(
        d_vector_size, profugus::Global_RNG::d_rng );

    // START PROFILING
    cudaProfilerStart();

    // run all the local histories while the source exists and there are live
    // particles in the vector. we know when all the particles are dead when
    // the starting point for dead events is at the front of the vector. there
    // is no need to communicate particles because the problem is replicated
    size_type dead_start = 0;
    size_type dead_end = 0;
    while ( !d_source->empty() || (0 != dead_start) )
    {
        // Run the events. Right now this works sequentially because there are
        // only 3 events - sampling the source, stepping to a collision or
        // boundary, and processing a collision or boundary. In the
        // future, we will first have to find the vector slots each event
        // kernel will operate on and then launch the kernels.
        d_source->get_particles( particles );
        d_transporter.transport_step( particles, bank );
        d_transporter.process_step( particles, bank );

        // Sort the vector.
        particles.get_host_ptr()->sort_by_event();

        // Find where the dead particles start.
        particles.get_host_ptr()->get_event_particles( 
            events::DEAD, dead_start, dead_end );

        // update the event counter
        ++counter;

        // print message if needed
        if (counter % d_print_count == 0)
        {
            double percent_complete
                = (100. * (d_source->num_run()-dead_start)) / (d_source->num_to_transport());
            cout << ">>> Finished " << counter << " events and ("
                 << std::setw(6) << std::fixed << std::setprecision(2)
                 << percent_complete << "%) particles on domain "
                 << d_node << endl;
        }
    }

    // barrier at the end
    profugus::global_barrier();

    // STOP PROFILING
    cudaProfilerStop();
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
    const SP_Fission_Sites& fis_sites,
    double keff)
{
    // set the transporter with the fission site container and the latest keff
    // iterate
    d_transporter.set(fis_sites, keff);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the variance reduction.
 */
template <class Geometry>
void Source_Transporter<Geometry>::set(SP_Variance_Reduction vr)
{
    REQUIRE(vr);

    // set the variance reduction in the domain transporter and locally
    d_transporter.set(vr);
    d_var_reduction = vr;

    ENSURE(d_var_reduction);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the tally controller
 */
template <class Geometry>
void Source_Transporter<Geometry>::set(const SP_Tallier& tallier)
{
    REQUIRE(tallier);

    // set the tally controller in the domain transporter and locally
    d_transporter.set(tallier);
    d_tallier = tallier;

    ENSURE(d_tallier);
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Source_Transporter_t_cuh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.t.hh
//---------------------------------------------------------------------------//
