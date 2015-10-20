//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Transporter.cc
 * \author Thomas M. Evans
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iomanip>
#include <iostream>
#include <cmath>

#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "Source_Transporter.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor/
 */
Source_Transporter::Source_Transporter(RCP_Std_DB  db,
                                       SP_Geometry geometry,
                                       SP_Physics  physics)
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
    d_print_fraction = db->get("mc_diag_frac", 1.1);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Assign the source.
 */
void Source_Transporter::assign_source(SP_Source source)
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
void Source_Transporter::solve()
{
    using std::cout; using std::endl;

    REQUIRE(d_source);

    // barrier at the start
    profugus::global_barrier();

    SCOPED_TIMER("MC::Source_Transporter.solve");

    // particle counter
    size_type counter = 0;

    // get a base class reference to the source
    Source_t &source = *d_source;

    // make a particle bank
    typename Transporter_t::Bank_t bank;
    CHECK(bank.empty());

    // run all the local histories while the source exists, there is no need
    // to communicate particles because the problem is replicated
    int np = source.num_to_transport();
    for ( int i = 0; i < np; ++i )
    {
        // get a particle from the source
        SP_Particle p = source.get_particle();
        CHECK(p);
        CHECK(p->alive());

        // Do "source event" tallies on the particle
        d_tallier->source(*p);

        // transport the particle through this (replicated) domain
        d_transporter.transport(*p, bank);
        CHECK(!p->alive());

        // transport any secondary particles that are part of this history
        // (from splitting or physics) that get put into the bank
        while (!bank.empty())
        {
            // get a particle from the bank
            SP_Particle bank_particle = bank.pop();
            CHECK(bank_particle);
            CHECK(bank_particle->alive());

            // make particle alive
            bank_particle->live();

            // transport it
            d_transporter.transport(*bank_particle, bank);
            CHECK(!bank_particle->alive());
        }

        // update the counter
        ++counter;

        // indicate completion of particle history
        d_tallier->end_history(*p);

        // print message if needed
        if (counter % d_print_count == 0)
        {
            double percent_complete
                = (100. * counter) / source.num_to_transport();
            cout << ">>> Finished " << counter << "("
                 << std::setw(6) << std::fixed << std::setprecision(2)
                 << percent_complete << "%) particles on domain "
                 << d_node << endl;
        }
    }

    // barrier at the end
    profugus::global_barrier();

    // check that we transported all of the particles.
    CHECK( source.empty() );

    // increment the particle counter
    DIAGNOSTICS_ONE(integers["particles_transported"] += counter);

#ifdef REMEMBER_ON
    profugus::global_sum(counter);
    ENSURE(counter == source.total_num_to_transport());
    ENSURE(bank.empty());
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the variance reduction.
 */
void Source_Transporter::set(SP_Variance_Reduction vr)
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
void Source_Transporter::set(SP_Tallier tallier)
{
    REQUIRE(tallier);

    // set the tally controller in the domain transporter and locally
    d_transporter.set(tallier);
    d_tallier = tallier;

    ENSURE(d_tallier);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.cc
//---------------------------------------------------------------------------//
