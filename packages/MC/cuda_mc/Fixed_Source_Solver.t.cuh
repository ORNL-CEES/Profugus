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
Fixed_Source_Solver<Geometry>::Fixed_Source_Solver(RCP_Std_DB db)
    : d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
    d_batch_size = std::numeric_limits<size_type>::max();
    if (db->isType<int>("batch_size"))
    {
        d_batch_size = db->get<int>("batch_size");
    }
    else if (db->isType<size_type>("batch_size"))
    {
        d_batch_size = db->get<size_type>("batch_size");
    }
    else if (db->isParameter("batch_size"))
        VALIDATE(false,"Unrecognized type for parameter batch_size.");
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the underlying source transporter and source.
 */
template <class Geometry>
void Fixed_Source_Solver<Geometry>::set(SP_Source_Transporter transporter,
                                        SP_Source             source,
                                        SP_Tallier_DMM        tallier)
{
    REQUIRE(transporter);
    REQUIRE(source);
    REQUIRE(tallier);

    // assign the solver
    d_transporter = transporter;

    // assign the source
    d_source = source;

    // assign tallier
    b_tallier = tallier;
    d_transporter->set(b_tallier);
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
    REQUIRE(b_tallier);

    SCOPED_TIMER("MC::Fixed_Source_Solver.solve");

    // store the number of source particles
    DIAGNOSTICS_TWO(vec_integers["local_np"].push_back(
                        d_source->num_to_transport()));

    // start the timer
    profugus::Timer fixed_source_timer;
    fixed_source_timer.start();

    // Determine source type
    auto uni_source =
        std::dynamic_pointer_cast<Uniform_Source_DMM<Geometry>>(d_source);
    if( uni_source )
    {
        d_Np = uni_source->total_num_to_transport();

        // Compute actual batch size
        uni_source->set_batch_size(d_batch_size);

        while (uni_source->num_left() > 0)
        {
            // Get actual particle count for this batch
            auto num_batch = uni_source->num_batch();

            // Transport batch of particles
            d_transporter->solve(uni_source);

            // Finalize batch in source and tallier
            b_tallier->end_batch(static_cast<double>(num_batch));
        }
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
    b_tallier->finalize(d_Np);
}

} // end namespace cuda_mc

#endif // cuda_mc_Fixed_Source_Solver_t_cuh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.t.cuh
//---------------------------------------------------------------------------//
