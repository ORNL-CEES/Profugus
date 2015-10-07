//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstDataFlow.cc
 * \author Stuart Slattery
 * \brief  HPX plain action testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>
#include <hpx/lcos/local/dataflow.hpp>

#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test Helpers
//---------------------------------------------------------------------------//
// Function for adding two numbers in the global namespace.
double add_numbers( const double a, const double b )
{
    return a + b;
}

// Function for multiplying two numbers in the global namespace.
double multiply_numbers( const double a, const double b )
{
    return a * b;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( plain_action, add_numbers_test )
{
    // Initialize values.
    double init_a = 1.0;
    double init_b = 2.0;
    hpx::shared_future<double> a = hpx::make_ready_future( init_a );
    hpx::shared_future<double> b = hpx::make_ready_future( init_b );

    // Run calculations. As soon as the futures passed to a function in a
    // dataflow are ready, the data flow function will execute automatically.
    int num_run = 100;
    for ( int i = 0; i < num_run; ++i )
    {
	hpx::shared_future<double> c = 
	    hpx::lcos::local::dataflow( 
		hpx::util::unwrapped(multiply_numbers), a, b );
	a = hpx::lcos::local::dataflow( 
	    hpx::util::unwrapped(add_numbers), a, c );
    }
    
    // Test the results.
    EXPECT_EQ( std::pow(init_a*(init_b+1),num_run), a.get() );
}

//---------------------------------------------------------------------------//
//                 end of tstDataFlow.cc
//---------------------------------------------------------------------------//
