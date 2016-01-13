//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstCollision_Tally.cc
 * \author Stuart Slattery
 * \brief  Collision_Tally class test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Collision_Tally.hh"
#include "Collision_Tally_Tester.hh"

#include "rng/RNG_Control.hh"

#include <Teuchos_Array.hpp>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(Collision_Tally, construction)
{
    // Geometry parameters.
    int N = 10;
    Teuchos::Array<double> x_edges( N+1 );
    Teuchos::Array<double> y_edges( N+1 );
    Teuchos::Array<double> z_edges( N+1 );
    for ( int i = 0; i < N+1; ++i )
    {
	x_edges[i] = i;
	y_edges[i] = i;
	z_edges[i] = i;
    }
    int num_cells = N*N*N;

    // Tally paramters.
    int num_batch = 30;

    // Particle vector parameters
    int num_particle = num_batch*num_cell;
    profugus::RNG_Control control( 3420239343 );

    // Create the tally tester.
    Collision_Tally_Tester tester( 
	x_edges, y_edges, z_edges, num_particle, control.rng(), num_batch );

    // Create a tally.
    Collision_Tally tally( tester.geometry(), num_batch );

    // check sizes
    EXPECT_EQ( tally.num_cells(), num_cells );
    EXPECT_EQ( tally.num_batch(), num_batch );

    // Tally the particles.
    tally.tally( tester.particles() );

    // Finalize the particles.
    tally.finalize( num_particle );

    // Get the moments.
    Teuchos::Array<double> first_moment;
    Teuchos::Array<double> second_moment;
    tally.copy_moments_to_host( first_moment, second_moment );
    EXPECT_EQ( first_moment.size(), num_cells );
    EXPECT_EQ( second_moment.size(), num_cells );

    // Check the first moment. Only the first half of the cells should have
    // particles with collisions that will tally.
    double gold_first = 0.0;
    for ( int b = 0; b < num_batch; ++b )
    {
	gold_first += i;
    }
    for ( int c = 0; c < num_cells; ++c )
    {
	if ( c < num_cells/2 )
	{
	    EXPECT_EQ( first_moment[c], c*gold_first/num_particle );
	}
	else
	{
	    EXPECT_EQ( first_moment[c], 0.0 );	    
	}
    }

    // Check the second moment.
    double gold_second = 0.0;
    for ( int b = 0; b < num_batch; ++b )
    {
	gold_second += i*i - 
		       gold_first*gold_first/(num_particle*num_particle);
    }
    for ( int c = 0; c < num_cells; ++c )
    {
	if ( c < num_cells/2 )
	{
	    EXPECT_EQ( second_moment[c], 
		       c*c*gold_second/(num_particle*(num_particle-1) );
	}
	else
	{
	    EXPECT_EQ( second_moment[c], 0.0 );	    
	}
    }
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
