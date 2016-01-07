//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstParticle.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 11:32:36 2014
 * \brief  Particle class test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Particle_Vector.hh"
#include "Particle_Vector_Tester.hh"

#include "mc/Definitions.hh"

#include "rng/RNG_Control.hh"

#include <algorithm>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(Particle, construction)
{
    int num_particle = 4;

    profugus::RNG_Control control( 3420239343 );

    Particle_Vector_Tester tester( num_particle, control.rng() );

    // check size
    EXPECT_EQ( tester.size(), num_particle );

    // check random number generation
    Teuchos::Array<double> random = tester.ran();
    for ( auto& r : random )
    {
	bool passed = ( (r >= 0.0) && (r <= 1.0) );
	EXPECT_TRUE( passed );
    }

    // check weight
    double wt = 1.34;
    tester.set_wt( wt );
    Teuchos::Array<double> weights = tester.wt();
    for ( auto& w : weights )
    {
	EXPECT_EQ( w, wt );
    }
    tester.multiply_wt( wt );
    weights = tester.wt();
    for ( auto& w : weights )
    {
	EXPECT_EQ( w, wt*wt );
    }

    // check group
    int grp = 32;
    tester.set_group( grp );
    Teuchos::Array<int> groups = tester.group();
    for ( auto& g : groups )
    {
	EXPECT_EQ( g, grp );
    }

    // Check matid.
    int mid = 19;
    tester.set_matid( mid );
    Teuchos::Array<int> matids = tester.matid();
    for ( auto& m : matids )
    {
	EXPECT_EQ( m, mid );
    }

    // Check alive status.
    tester.live();
    Teuchos::Array<int> alive = tester.alive();
    for ( auto& a : alive )
    {
	EXPECT_TRUE( a );
    }
    tester.kill();
    alive = tester.alive();
    for ( auto& a : alive )
    {
	EXPECT_FALSE( a );
    }

    // Setup events.
    typedef typename Particle_Vector_Tester::Event_t Event_t;
    Teuchos::Array<Event_t> host_events( num_particle );
    for ( int i = 0; i < num_particle; ++i )
    {
	// Evens scatter, odds absorb
	host_events[i] = ( i % 2 == 0 ) 
			 ? profugus::events::SCATTER 
			 : profugus::events::ABSORPTION;
    }

    // Check event assignment
    tester.set_event( host_events );
    Teuchos::Array<Event_t> device_events = tester.event();
    for ( int i = 0; i < num_particle; ++i )
    {
	EXPECT_EQ( host_events[i], device_events[i] );
    }

    // Check event sorting.
    tester.sort_by_event();
    std::sort( host_events.begin(), host_events.end() );
    device_events = tester.event();
    for ( int i = 0; i < num_particle; ++i )
    {
	EXPECT_EQ( host_events[i], device_events[i] );
    }

    // Setup a geo state.
    typedef typename Particle_Vector_Tester::Geo_State_t Geo_State_t;
    Geo_State_t geo_state;

    geo_state.ijk.i = 2.3;
    geo_state.ijk.j = 1.3;
    geo_state.ijk.k = 3.3;

    geo_state.d_r.x = 2.4;
    geo_state.d_r.y = 1.4;
    geo_state.d_r.z = 3.4;

    geo_state.d_dir.x = 2.5;
    geo_state.d_dir.y = 1.5;
    geo_state.d_dir.z = 3.5;

    geo_state.next_ijk.i = 2.6;
    geo_state.next_ijk.j = 1.6;
    geo_state.next_ijk.k = 3.6;

    geo_state.next_dist = 4.3;

    // check the geo state
    tester.set_geo_state( geo_state );
    Teuchos::Array<Geo_State_t> states = tester.geo_state();
    for ( auto& s : states )
    {
	EXPECT_EQ( geo_state.ijk.i, s.ijk.i );
	EXPECT_EQ( geo_state.ijk.j, s.ijk.j );
	EXPECT_EQ( geo_state.ijk.k, s.ijk.k );

	EXPECT_EQ( geo_state.d_r.x, s.d_r.x );
	EXPECT_EQ( geo_state.d_r.y, s.d_r.y );
	EXPECT_EQ( geo_state.d_r.z, s.d_r.z );

	EXPECT_EQ( geo_state.d_dir.x, s.d_dir.x );
	EXPECT_EQ( geo_state.d_dir.y, s.d_dir.y );
	EXPECT_EQ( geo_state.d_dir.z, s.d_dir.z );

	EXPECT_EQ( geo_state.next_ijk.i, s.next_ijk.i );
	EXPECT_EQ( geo_state.next_ijk.j, s.next_ijk.j );
	EXPECT_EQ( geo_state.next_ijk.k, s.next_ijk.k );

	EXPECT_EQ( geo_state.next_dist, s.next_dist );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
