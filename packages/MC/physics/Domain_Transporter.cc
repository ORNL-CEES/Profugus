//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Domain_Transporter.cc
 * \author Thomas M. Evans
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "harness/DBC.hh"
#include "harness/Diagnostics.hh"
#include "geometry/Definitions.hh"
#include "Definitions.hh"
#include "Domain_Transporter.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Domain_Transporter::Domain_Transporter()
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the geometry and physics classes.
 *
 * \param geometry
 * \param physics
 */
void Domain_Transporter::set(SP_Geometry geometry,
                             SP_Physics  physics)
{
    REQUIRE(geometry);
    REQUIRE(physics);

    d_geometry = geometry;
    d_physics  = physics;

    if (d_var_reduction)
    {
        d_var_reduction->set(d_geometry);
        d_var_reduction->set(d_physics);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the variance reduction.
 *
 * \param reduction
 */
void Domain_Transporter::set(SP_Variance_Reduction reduction)
{
    REQUIRE(reduction);
    d_var_reduction = reduction;

    if (d_geometry)
        d_var_reduction->set(d_geometry);
    if (d_physics)
        d_var_reduction->set(d_physics);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set regular tallies.
 *
 * \param tallies
 */
void Domain_Transporter::set(SP_Tallier tallies)
{
    REQUIRE(tallies);
    d_tallier = tallies;
    ENSURE(d_tallier);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set regular source.
 *
 * \param source
 */
void Domain_Transporter::set(SP_Source source)
{
    REQUIRE(source);
    d_source = source;
    ENSURE(d_source);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transport a particle through the domain.
 *
 * \param particle
 * \param bank
 */
void Domain_Transporter::transport( 
    std::vector<Particle_t>& particles,
    std::vector<std::pair<std::size_t,events::Event> >& events, 
    std::vector<Bank_t>& banks ) const
{
    REQUIRE(d_geometry);
    REQUIRE(d_physics);
    REQUIRE(d_var_reduction);
    REQUIRE(d_tallier);

    // Get the number of particles
    int batch_size = particles.size();
    CHECK( batch_size > 0 );

    // initialize the source run count.
    def::size_type num_run = batch_size;
    def::size_type num_to_run = d_source->num_to_transport();

    // initialize the particle distances.
    for ( int n = 0; n < batch_size; ++n )
    {
	particles[n].set_dist_mfp( -std::log(particles[n].rng().ran()) );
    }

    // Create an event sorting and searching predicate.
    auto event_pred = []( const std::pair<std::size_t,events::Event>& e1,
			  const std::pair<std::size_t,events::Event>& e2 )
		      { return e1.second < e2.second; };

    // Create functions for find the range of collision events.
    auto find_collision_begin = []( const std::pair<std::size_t,events::Event>& e )
				{ return e.second >= events::COLLISION; };
    auto find_collision_end = []( const std::pair<std::size_t,events::Event>& e )
				{ return e.second > events::COLLISION; };

    // Create functions for find the range of boundary events.
    auto find_boundary_begin = []( const std::pair<std::size_t,events::Event>& e )
				{ return e.second >= events::BOUNDARY; };
    auto find_boundary_end = []( const std::pair<std::size_t,events::Event>& e )
				{ return e.second > events::BOUNDARY; };

    // Create functions for find the range of end history events.
    auto find_end_history_begin = []( const std::pair<std::size_t,events::Event>& e )
				{ return e.second > events::STILL_ALIVE; };
    auto find_end_history_end = []( const std::pair<std::size_t,events::Event>& e )
				{ return e.second == events::END_EVENT; };
    
    // Create an iterator to the end of the active particles. The whole array
    // is alive to start.
    std::vector<std::pair<std::size_t,events::Event> >::iterator alive_end
	= events.end();
    std::vector<std::pair<std::size_t,events::Event> >::iterator dead_end
	= events.end();

    // Run the kernels until the first particle in the list is no longer
    // alive. The sorting will make sure that when this is the case, all
    // particles have been finished.
    int counter = 0;
    while ( events[0].second < events::STILL_ALIVE )
    {
	// Get the next event for the particles that are still alive.
	for ( auto e = events.begin(); e != alive_end; ++e )
	{
	    get_next_event( particles[e->first], e->second );
	}

	// Sort the alive particle events.
	std::sort( events.begin(), alive_end, event_pred );

	// Get the range of particles that have had a collision.
	auto collision_range_begin = std::find_if( events.begin(),
						   alive_end,
						   find_collision_begin );
	auto collision_range_end = std::find_if( events.begin(),
						 alive_end,
						 find_collision_end );

	// Get the range of particles that have hit a boundary.
	auto boundary_range_begin = std::find_if( events.begin(),
						  alive_end,
						  find_boundary_begin );
	auto boundary_range_end = std::find_if( events.begin(),
						alive_end,
						find_boundary_end );

	// Process collision events.
	for ( auto e = collision_range_begin; e != collision_range_end; ++e )
	{
	    process_collision( particles[e->first], e->second, banks[e->first] );
	    particles[e->first].set_dist_mfp(
		-std::log(particles[e->first].rng().ran()) );
	}

	// Process boundary events.
	for ( auto e = boundary_range_begin; e != boundary_range_end; ++e )
	{
	    process_boundary( particles[e->first], e->second, banks[e->first] );
	}

	// Tally particles that died on the last event cycle and spawn a new
	// one if needed.
	for ( auto e = alive_end; e != dead_end; ++e )
	{
	    d_tallier->end_history( particles[e->first] );

	    // If we still have particles to run create a
	    // new one.
	    if ( num_run < num_to_run )
	    {
		particles[e->first] = d_source->get_particle(num_run);
		e->second = events::BORN;
		++num_run;
	    }

	    // Otherwise this element in the vector is done.
	    else
	    {
		e->second = events::END_EVENT;
	    }
	}

	// Sort the events again to get the active particles and dead particles.
	std::sort( events.begin(), dead_end, event_pred );

	// Get the iterator to the end of the particles that are alive.
	alive_end = std::find_if( events.begin(),
				  dead_end,
				  find_end_history_begin );

	// Get the iterator to the beginning of the particles that are
	// completely dead.
	dead_end = std::find_if( events.begin(),
				 dead_end,
				 find_end_history_end );

	// Increment cycle counter.
	++counter;
    } 

    CHECK( alive_end == events.begin() );

    // Process any remaining particles.
    for ( auto e = alive_end; e != dead_end; ++e )
    {
	d_tallier->end_history( particles[e->first] );
    } 
    
    std::cout << "NUM RUN " << num_run << std::endl;
    std::cout << "NUM CYCLE " << counter << std::endl;
    ENSURE( num_run == num_to_run );
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
void Domain_Transporter::get_next_event( 
    Particle_t& particle, events::Event& event ) const
{
    Step_Selector step_selector;

    CHECK(particle.alive());
    CHECK(particle.dist_mfp() > 0.0);

    // total interaction cross section
    double xs_total = d_physics->total(physics::TOTAL, particle);
    CHECK(xs_total >= 0.0);

    // sample distance to next collision
    double dist_col = 
	(xs_total > 0.0) ? (particle.dist_mfp() / xs_total) : constants::huge;

    // initialize the distance to collision in the step selector
    step_selector.initialize(dist_col, events::COLLISION);

    // calculate distance to next geometry boundary
    double dist_bnd = d_geometry->distance_to_boundary( particle.geo_state() );
    step_selector.submit(dist_bnd, events::BOUNDARY);

    // set the next event in the particle
    CHECK(step_selector.tag() < events::END_EVENT);
    event = static_cast<events::Event>( step_selector.tag() );

    // path-length tallies (the actual movement of the particle will
    // take place when we process the various events)
    d_tallier->path_length(step_selector.step(), particle);

    // update the mfp distance travelled
    if (events::COLLISION == event)
    {
	particle.set_dist_mfp( step_selector.step() );
    }
    else
    {
	particle.add_dist_mfp( -step_selector.step()*xs_total );
    }
}

//---------------------------------------------------------------------------//
void Domain_Transporter::process_boundary(Particle_t &particle,
					  events::Event& event,
                                          Bank_t     &bank) const
{
    REQUIRE(particle.alive());
    REQUIRE(event == events::BOUNDARY);

    // return if not a boundary event

    // move a particle to the surface and process the event through the surface
    d_geometry->move_to_surface(particle.geo_state());

    // get the in/out state of the particle
    int state = d_geometry->boundary_state(particle.geo_state());

    // reflected flag
    bool reflected = false;

    // process the boundary crossing based on the geometric state of the
    // particle
    switch (state)
    {
        case geometry::OUTSIDE:
            // the particle has left the problem geometry
            event = events::ESCAPE;
            particle.kill();

            // add a escape diagnostic
            DIAGNOSTICS_TWO(integers["geo_escape"]++);

            ENSURE(event == events::ESCAPE);
            break;

        case geometry::REFLECT:
            // the particle has hit a reflecting surface
            reflected = d_geometry->reflect(particle.geo_state());
            CHECK(reflected);

            // add a reflecting face diagnostic
            DIAGNOSTICS_TWO(integers["geo_reflect"]++);

            ENSURE(event == events::BOUNDARY);
            break;

        case geometry::INSIDE:
            // otherwise the particle is at an internal geometry boundary;
            // update the material id of the region the particle has entered
            particle.set_matid(d_geometry->matid(particle.geo_state()));

            // add variance reduction at surface crossings
            d_var_reduction->post_surface(particle, event, bank);

            // add a boundary crossing diagnostic
            DIAGNOSTICS_TWO(integers["geo_surface"]++);

            ENSURE(event == events::BOUNDARY);
            break;

        default:
            CHECK(0);
    }
}

//---------------------------------------------------------------------------//

void Domain_Transporter::process_collision(Particle_t &particle,
					   events::Event& event,
                                           Bank_t     &bank) const
{
    REQUIRE(d_var_reduction);
    REQUIRE(event == events::COLLISION);

    // move the particle to the collision site
    d_geometry->move_to_point(particle.dist_mfp(), particle.geo_state());

    // use the physics package to process the collision
    d_physics->collide(particle, event, bank);

    // apply weight windows
    d_var_reduction->post_collision(particle, event, bank);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.cc
//---------------------------------------------------------------------------//
