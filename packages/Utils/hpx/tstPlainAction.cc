//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstPlainAction.cc
 * \author Stuart Slattery
 * \brief  HPX plain action testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>

#include <mutex>
#include <vector>
#include <algorithm>
#include <string>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test Helpers
//---------------------------------------------------------------------------//
// Function for adding two numbers in the global namespace.
double add_numbers( const double a, const double b )
{
    return a + b;
}

// Make a plain action out of the function. This is effectively a
// packaged_task but more nuanced.
HPX_PLAIN_ACTION( add_numbers, add_numbers_action );

//---------------------------------------------------------------------------//
// Future-based add function.
double add_futures( hpx::future<double> a, hpx::future<double> b )
{
    return a.get() + b.get();
}

HPX_PLAIN_ACTION( add_futures, add_futures_action );

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( plain_action, add_numbers_test )
{
    // Create test data.
    double a = 3.43;
    double b = -2.39999;

    // Create an instance of the action.
    add_numbers_action ana;

    // Execute the action on the local machine.
    double result = ana( hpx::find_here(), a, b );

    // Test the results.
    EXPECT_EQ( a+b, result );
}

//---------------------------------------------------------------------------//
TEST( plain_action, add_futures_test )
{
    // Create test data.
    double a = 3.43;
    double b = -2.39999;
    double c = 5.2;
    double d = 1.993;

    // Create an instance of the actions.
    add_numbers_action ana;
    add_futures_action afa;

    // Execute the action on the local machine.
    hpx::future<double> f1 = hpx::async( ana, hpx::find_here(), a, b );
    hpx::future<double> f2 = hpx::async( ana, hpx::find_here(), c, d );
    hpx::future<double> result = 
	hpx::async( afa, hpx::find_here(), std::move(f1), std::move(f2) );

    // Test the results.
    EXPECT_EQ( a+b+c+d, result.get() );
}

//---------------------------------------------------------------------------//
TEST( plain_action, wait_all_test )
{
    // Allocate space for futures of the local sum tasks.
    int num_data = 10000;
    std::vector<hpx::future<double> > local_sum_futures;
    local_sum_futures.reserve( num_data );

    // Initialize the tasks.
    for ( int i = 0; i < num_data; ++i )
    {
	local_sum_futures.push_back( 
	    hpx::async<add_numbers_action>( hpx::find_here(), i, i )
	    );
    }

    // Wait on the tasks, extract their data, and write into the global sum
    // with a mutex protection.
    hpx::lcos::local::spinlock mtx;
    double global_sum = 0.0;
    auto sum_func = [&](double s)
		    { 
			std::lock_guard<decltype(mtx)> lock(mtx);
			global_sum += s; 
		    };
    hpx::lcos::wait_each( hpx::util::unwrapped(sum_func), local_sum_futures );

    // Test the result.
    double gold_sum = 0.0;
    for ( int i = 0; i < num_data; ++i )
    {
	gold_sum += i + i;
    }
    EXPECT_EQ( global_sum, gold_sum );
}

//---------------------------------------------------------------------------//
//                 end of tstPlainAction.cc
//---------------------------------------------------------------------------//
