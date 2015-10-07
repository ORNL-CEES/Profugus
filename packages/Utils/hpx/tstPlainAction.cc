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
#include <hpx/include/util.hpp>

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
//                 end of tstPlainAction.cc
//---------------------------------------------------------------------------//
