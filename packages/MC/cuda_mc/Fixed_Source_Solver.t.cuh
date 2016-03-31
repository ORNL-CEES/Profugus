//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fixed_Source_Solver.t.cuh
 * \author Steven Hamilton
 * \date   Tue May 13 14:40:06 2014
 * \brief  Fixed_Source_Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fixed_Source_Solver_t_cuh
#define cuda_mc_Fixed_Source_Solver_t_cuh

#include "Fixed_Source_Solver.hh"

#include "Uniform_Source.cuh"
#include "Tallier.cuh"

#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/P_Stream.hh"
#include "comm/Timing.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Fixed_Source_Solver<Geometry>::Fixed_Source_Solver()
    : d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the underlying source transporter and source.
 */
template <class Geometry>
void Fixed_Source_Solver<Geometry>::set(SP_Source_Transporter transporter,
                                        SP_Source             source)
{
    REQUIRE(transporter);
    REQUIRE(source);

    // assign the solver
    d_transporter = transporter;

    // assign the source
    d_source = source;

    // get the tallies and assign them
    b_tallier = d_transporter->tallier();
    INSIST(b_tallier.get_host_ptr(),
            "Tally not assigned in Source_Transporter in fixed-source solver.");
    INSIST(b_tallier.get_device_ptr(),
            "Tally not assigned in Source_Transporter in fixed-source solver.");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 *
 * This also finalizes tallies.
 */
template <class Geometry>
void Fixed_Source_Solver<Geometry>::solve()
{
    using profugus::endl; using profugus::pcout;

    REQUIRE(d_source);
    REQUIRE(b_tallier.get_host_ptr());
    REQUIRE(b_tallier.get_device_ptr());

    SCOPED_TIMER("MC::Fixed_Source_Solver.solve");

    // store the number of source particles
    DIAGNOSTICS_TWO(vec_integers["local_np"].push_back(
                        d_source->num_to_transport()));

    // start the timer
    profugus::Timer fixed_source_timer;
    fixed_source_timer.start();

    // Determine source type
    auto uni_source =
        std::dynamic_pointer_cast<Uniform_Source<Geometry>>(d_source);
    if( uni_source )
    {
        d_Np = uni_source->total_num_to_transport();
        d_transporter->solve(uni_source);
    }

    // Add other source types...
    REQUIRE(uni_source);

    // stop the timer
    profugus::global_barrier();
    fixed_source_timer.stop();

    // output
    pcout << ">>> Finished transporting " << profugus::setw(8)
          << d_Np  << " particles in "
          << profugus::fixed << profugus::setw(12)
          << profugus::setprecision(6) << fixed_source_timer.TIMER_CLOCK()
          << " seconds" << endl;

    // Finalize tallies using global number of particles
    b_tallier.get_host_ptr()->finalize(d_Np);
}

} // end namespace cuda_mc

#endif // cuda_mc_Fixed_Source_Solver_t_cuh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.t.cuh
//---------------------------------------------------------------------------//