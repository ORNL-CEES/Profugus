//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstPlainAction.cc
 * \author Stuart Slattery
 * \brief  HPX plain action testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <string>

#include "gtest/utils_gtest_hpx.hh"

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
// TESTS
//---------------------------------------------------------------------------//
TEST( plain_action, add_test_1 )
{
    // Create 2 numbers.
    double a = 3.43;
    double b = -2.39999;

    // Create an instance of the action.
    add_numbers_action add;

    // Execute the action on the local machine.
    double result = add( hpx::find_here(), a, b );

    // Test the results.
    EXPECT_EQ( a+b, result );
}

//---------------------------------------------------------------------------//
//                 end of tstPlainAction.cc
//---------------------------------------------------------------------------//
