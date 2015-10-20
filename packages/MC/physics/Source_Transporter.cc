//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Transporter.cc
 * \author Thomas M. Evans
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <hpx/parallel/execution_policy.hpp>
#include <hpx/include/parallel_algorithm.hpp>

#include <boost/range/irange.hpp>

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
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 */
void Source_Transporter::solve()
{
    using std::cout; using std::endl;

    REQUIRE(d_source);

    SCOPED_TIMER("MC::Source_Transporter.solve");

    // run all the local histories while the source exists, there is no need
    // to communicate particles because the problem is replicated
    int np = d_source->num_to_transport();
    auto range = boost::irange( 0, np );
    auto transport_func = [this](const int){ this->transport_history(); };
    auto transport_task =
	hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 transport_func );
    transport_task.wait();

    // check that we transported all of the particles.
    CHECK( d_source->empty() );
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

//---------------------------------------------------------------------------//
/*! 
 * \brief Transport a history in the source and any histories it may make from
 * splitting.
 */
void Source_Transporter::transport_history()
{
    // make a particle bank
    typename Transporter_t::Bank_t bank;
    CHECK(bank.empty());

    // get a particle from the source
    SP_Particle p = d_source->get_particle();
    CHECK(p);
    CHECK(p->alive());

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

    // indicate completion of particle history
    d_tallier->end_history(*p);

    // Make sure the bank is empty.
    ENSURE(bank.empty());
}

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.cc
//---------------------------------------------------------------------------//
