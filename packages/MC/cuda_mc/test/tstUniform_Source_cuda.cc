//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstUniform_Source_cuda.cc
 * \author Stuart Slattery
 * \brief  Uniform_Source class test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Uniform_Source.hh"
#include "../Box_Shape.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "geometry/Definitions.hh"

#include "Uniform_Source_Tester.hh"

#include "rng/RNG_Control.hh"

#include <Teuchos_Array.hpp>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(Uniform_Source, construction)
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
    typedef typename profugus::geometry::matid_type matid_type;
    matid_type matid = 3;

    // Particle vector parameters
    int vector_size = 256;
    profugus::RNG_Control control( 3420239343 );

    // Source parameters.
    int num_batch = 2;
    int num_source = num_batch*vector_size;
    int num_group = 4;
    Teuchos::Array<double> spectral_shape( num_group, 1.0 );
    Teuchos::RCP<Teuchos::ParameterList> db = Teuchos::parameterList();
    db->set( "Np", num_source );
    db->set( "spectral_shape", spectral_shape );

    // Create the source tester.
    Uniform_Source_Tester tester( x_edges, y_edges, z_edges, matid, 
				  vector_size, control.rng(),
				  num_group );

    // Create a source.
    cuda_profugus::Uniform_Source<cuda_profugus::Mesh_Geometry,
				  cuda_profugus::Box_Shape> 
	source( db, tester.geometry(), num_group, num_batch );
    
    // Build the source.
    source.build_source( tester.shape() );

    // Check the source accessors prior to transport.
    EXPECT_FALSE( source.empty() );
    EXPECT_EQ( source.num_batch(), num_batch );
    EXPECT_EQ( source.num_to_transport(), num_source / profugus::nodes() );
    EXPECT_EQ( source.total_num_to_transport(), num_source );
    EXPECT_EQ( source.Np(), num_source );
    EXPECT_EQ( source.num_run(), 0 );
    EXPECT_EQ( source.num_left(), num_source / profugus::nodes() );

    // Run the first batch of particles.
    auto particles = tester.particles();
    source.get_particles( particles );

    // Check the state of the particles.
    Teuchos::Array<int> matids = tester.matid();
    Teuchos::Array<double> wts = tester.wt();
    Teuchos::Array<int> alive = tester.alive();
    Teuchos::Array<int> groups = tester.group();
    Teuchos::Array<cuda_profugus::events::Event> events = tester.event();
    Teuchos::Array<int> batches = tester.batch();
    for ( int i = 0; i < vector_size; ++i )
    {
    	EXPECT_EQ( matids[i], matid );
    	EXPECT_EQ( wts[i], 1.0 );
    	EXPECT_TRUE( alive[i] );
    	EXPECT_TRUE( (groups[i] >= 0) && (groups[i] < num_group) );
    	EXPECT_EQ( events[i], cuda_profugus::events::BORN );
    	EXPECT_EQ( batches[i], 0 );
    }

    // Check the state of the source.
    EXPECT_FALSE( source.empty() );
    EXPECT_EQ( source.num_run(), num_source / (2*profugus::nodes()) );
    EXPECT_EQ( source.num_left(), num_source / (2*profugus::nodes()) );

    // Kill all the particles so we can get a new batch.
    tester.kill_particles();

    // Run the second batch.
    source.get_particles( particles );

    // Check the state of the particles.
    matids = tester.matid();
    wts = tester.wt();
    alive = tester.alive();
    groups = tester.group();
    events = tester.event();
    batches = tester.batch();
    for ( int i = 0; i < vector_size; ++i )
    {
    	EXPECT_EQ( matids[i], matid );
    	EXPECT_EQ( wts[i], 1.0 );
    	EXPECT_TRUE( alive[i] );
    	EXPECT_TRUE( (groups[i] >= 0) && (groups[i] < num_group) );
    	EXPECT_EQ( events[i], cuda_profugus::events::BORN );
    	EXPECT_EQ( batches[i], 1 );
    }

    // Check that the source is empty.
    EXPECT_TRUE( source.empty() );
    EXPECT_EQ( source.num_run(), num_source / profugus::nodes() );
    EXPECT_EQ( source.num_left(), 0 );
}

//---------------------------------------------------------------------------//
//                 end of tstUniform_Source_cuda.cc
//---------------------------------------------------------------------------//
