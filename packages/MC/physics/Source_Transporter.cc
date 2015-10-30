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
    d_transporter.set( d_source );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the fixed-source problem.
 */
void Source_Transporter::solve()
{
    REQUIRE(d_source);

    SCOPED_TIMER("MC::Source_Transporter.solve");

    // Build a vector of particles.
    int batch_size = d_source->batch_size();
    CHECK( batch_size > 0 );
    std::vector<Particle_t> particles( batch_size );
    auto range = boost::irange( 0, batch_size );
    auto source_func = [&,this]( const int lid )
		       { particles[lid] = this->d_source->get_particle(lid); };
    auto source_task = hpx::parallel::for_each( 
	hpx::parallel::parallel_task_execution_policy(),
	std::begin(range),
	std::end(range),
	source_func );

    // Make a vector of event references.
    std::vector<std::pair<std::size_t,events::Event> > events( batch_size );
    auto event_fill_func = [&]( const int lid ) 
			   { events[lid] = std::make_pair(lid,events::BORN); };
    auto event_fill_task = hpx::parallel::for_each( 
	hpx::parallel::parallel_task_execution_policy(),
	std::begin(range),
	std::end(range),
	event_fill_func );

    // Make a vector of banks.
    std::vector<Bank_t> banks( batch_size );    

    // Wait to finish making particles and events.
    source_task.wait();
    event_fill_task.wait();

    // Transport all particles in the vector.
    d_transporter.transport( particles, events, banks );

    // Check that all of the banks are empty.
    REMEMBER( for ( auto& b : banks ) CHECK( b.empty() ); );
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
 * \brief Set the tallier.
 */
void Source_Transporter::set(SP_Tallier tallier)
{
    REQUIRE(tallier);

    // set the tally in the domain transporter and locally
    d_transporter.set(tallier);
    d_tallier = tallier;

    ENSURE(d_tallier);
}

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.cc
//---------------------------------------------------------------------------//
