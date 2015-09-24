//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstLocal_Smart_Thread.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 9:49:03 2015
 * \brief  Local_Smart_Thread testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Local_Smart_Thread.hh"

#include <vector>
#include <algorithm>
#include <thread>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Local_Smart_Thread, function1)
{
    // Declare test data outside of local thread scope.
    int first_odd = -1;
    int first_even = -1;

    // Calculate the test data.
    {
	std::vector<int> x = {10, 22, 31, 44, 56};

	// run a thread to find first odd
	auto odd_func = []( int& i, std::vector<int>& vec )
			{
			    i =
			    *std::find_if(std::begin(vec), std::end(vec),
					  [](int n){ return n % 2 == 1; });
			};
	profugus::Local_Smart_Thread odd_thread(
	    std::thread(odd_func, std::ref(first_odd), std::ref(x)) );

	// run a thread to find first even
	auto even_func = []( int& i, std::vector<int>& vec )
			 {
			     i =
			     *std::find_if(std::begin(vec), std::end(vec),
					   [](int n){ return n % 2 == 0; });
			 };

	profugus::Local_Smart_Thread even_thread =
	    profugus::make_lst(
		std::thread(even_func, std::ref(first_even), std::ref(x)) );
    
	// check that the threads have different ids
	EXPECT_NE( odd_thread.get_id(), even_thread.get_id() );
    }
    
    // check the even and odd results after the threads have been joined by
    // exiting the local scope to ensure the computation has been completed.
    EXPECT_EQ(31, first_odd);
    EXPECT_EQ(10, first_even );
}

//---------------------------------------------------------------------------//
//                 end of tstLocal_Smart_Thread.cc
//---------------------------------------------------------------------------//
