//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstCollision_Tally_cuda.cc
 * \author Stuart Slattery
 * \brief  Collision_Tally class test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Collision_Tally.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

#include "Tally_Tester.hh"

#include "rng/RNG_Control.hh"

#include <Teuchos_Array.hpp>
#include <Teuchos_TwoDArray.hpp>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(Collision_Tally, construction)
{
    // Geometry parameters.
    int N = 2;
    std::vector<double> x_edges( N+1 );
    std::vector<double> y_edges( N+1 );
    std::vector<double> z_edges( N+1 );
    for ( int i = 0; i < N+1; ++i )
    {
	x_edges[i] = i;
	y_edges[i] = i;
	z_edges[i] = i;
    }
    int num_cells = N*N*N;

    // Tally parameters.
    int num_batch = 256;

    // Particle vector parameters
    int num_particle = num_batch*num_cells;
    profugus::RNG_Control control( 3420239343 );

    // Create the tally tester.
    Tally_Tester tester( 
	x_edges, y_edges, z_edges, num_particle, control.rng(), num_batch );

    // Create a tally.
    cuda_profugus::Collision_Tally<cuda_profugus::Mesh_Geometry> 
	tally( tester.geometry(), num_batch );

    // check sizes
    EXPECT_EQ( tally.num_cells(), num_cells );
    EXPECT_EQ( tally.num_batch(), num_batch );

    // Tally the particles.
    tally.accumulate( tester.particles() );

    // Finalize the particles.
    tally.finalize( num_particle );

    // Get the moments.
    Teuchos::Array<double> first_moment;
    Teuchos::Array<double> second_moment;
    tally.copy_moments_to_host( first_moment, second_moment );
    EXPECT_EQ( first_moment.size(), num_cells );
    EXPECT_EQ( second_moment.size(), num_cells );

    // Calcuate the expected tally results. Only the first half of the cells
    // should have particles with collisions that will tally.
    Teuchos::TwoDArray<double> gold_tally( num_batch, num_cells, 0.0 );
    for ( int c = 0; c < num_cells; ++c )
    {
	if ( c < num_cells / 2 )
	{
	    for ( int b = 0; b < num_batch; ++b )
	    {
		gold_tally(b,c) += (b+1)*(c+1);
	    }
	}
    }
    for ( int c = 0; c < num_cells; ++c )
    {
	for ( int b = 0; b < num_batch; ++b )
	{
	    gold_tally(b,c) = gold_tally(b,c) * num_batch / num_particle;
	}
    }

    Teuchos::Array<double> gold_first( num_cells, 0.0 );
    for ( int c = 0; c < num_cells; ++c )
    {
	for ( int b = 0; b < num_batch; ++b )
	{
	    gold_first[c] += gold_tally(b,c) / num_batch;
	}
    }

    Teuchos::Array<double> gold_second( num_cells, 0.0 );
    for ( int c = 0; c < num_cells; ++c )
    {
	for ( int b = 0; b < num_batch; ++b )
	{
	    gold_second[c] += (gold_tally(b,c) * gold_tally(b,c) -
			       gold_first[c] * gold_first[c]) /
			      ( num_batch * (num_batch-1) );
	}
    }

    // Check that the tallies are correct.
    for ( int c = 0; c < num_cells; ++c )
    {
	EXPECT_FLOAT_EQ( gold_first[c], first_moment[c] );
	EXPECT_FLOAT_EQ( gold_second[c], second_moment[c] );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstCollision_Tally_cuda.cc
//---------------------------------------------------------------------------//
