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
    Require (!db.is_null());
    Require (d_geometry);
    Require (d_physics);

    // set the geometry and physics in the domain transporter
    d_transporter.set(d_geometry, d_physics);

    // set the output frequency for particle transport diagnostics
    d_print_fraction = db->get("mc_diag_frac", 1.0);
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

    Require (source);

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

    Require (d_source);

    // barrier at the start
    profugus::global_barrier();

    SCOPED_TIMER("MC::Source_Transporter.solve");

    // particle counter
    size_type counter = 0;

    // get a base class reference to the source
    Source_t &source = *d_source;

    // make a particle bank
    typename Transporter_t::Bank_t bank;
    Check (bank.empty());

    // run all the local histories while the source exists, there is no need
    // to communicate particles because the problem is replicated
    while (!source.empty())
    {
        // get a particle from the source
        SP_Particle p = source.get_particle();
        Check (p);
        Check (p->alive());

        // Do "source event" tallies on the particle
        d_tallier->source(*p);

        // transport the particle through this (replicated) domain
        d_transporter.transport(*p, bank);
        Check (!p->alive());

        // transport any secondary particles that are part of this history
        // (from splitting or physics) that get put into the bank
        while (!bank.empty())
        {
            // get a particle from the bank
            SP_Particle bank_particle = bank.pop();
            Check (bank_particle);
            Check (bank_particle->alive());

            // make particle alive
            bank_particle->live();

            // transport it
            d_transporter.transport(*bank_particle, bank);
            Check (!bank_particle->alive());
        }

        // update the counter
        ++counter;

        // indicate completion of particle history
        d_tallier->end_history();

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

    // increment the particle counter
    DIAGNOSTICS_ONE(integers["particles_transported"] += counter);

#ifdef REMEMBER_ON
    profugus::global_sum(counter);
    Ensure (counter == source.total_num_to_transport());
    Ensure (bank.empty());
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
void Source_Transporter::sample_fission_sites(SP_Fission_Sites fis_sites,
                                              double           keff)
{
    // set the transporter with the fission site container and the latest keff
    // iterate
    d_transporter.set(fis_sites, keff);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the variance reduction.
 */
void Source_Transporter::set(SP_Variance_Reduction vr)
{
    Require (vr);

    // set the variance reduction in the domain transporter and locally
    d_transporter.set(vr);
    d_var_reduction = vr;

    Ensure (d_var_reduction);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the tally controller
 */
void Source_Transporter::set(SP_Tallier tallier)
{
    Require (tallier);

    // set the tally controller in the domain transporter and locally
    d_transporter.set(tallier);
    d_tallier = tallier;

    Ensure (d_tallier);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.cc
//---------------------------------------------------------------------------//
