//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstThread.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 9:49:03 2015
 * \brief  Thread testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <thread>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(thread, thread_test)
{
    
    std::vector<int> x = {10, 22, 31, 44, 56};

    // run a thread to find first odd
    int first_odd = -1;
    auto odd_func = []( int& i, std::vector<int>& vec )
		    {
			i =
			*std::find_if(std::begin(vec), std::end(vec),
				      [](int n){ return n % 2 == 1; });
		    };
    std::thread odd_thread( odd_func, std::ref(first_odd), std::ref(x) );

    // run a thread to find first even
    int first_even = -1;
    auto even_func = []( int& i, std::vector<int>& vec )
		     {
			 i =
			 *std::find_if(std::begin(vec), std::end(vec),
				       [](int n){ return n % 2 == 0; });
		     };
    std::thread even_thread( even_func, std::ref(first_even), std::ref(x) );
    
    // check that the threads have different ids
    EXPECT_NE( odd_thread.get_id(), even_thread.get_id() );
    
    // join the threads before exiting
    EXPECT_TRUE( odd_thread.joinable() );
    EXPECT_TRUE( even_thread.joinable() );
    odd_thread.join();
    even_thread.join();
    EXPECT_FALSE( odd_thread.joinable() );
    EXPECT_FALSE( even_thread.joinable() );

    // check the even and odd results after the threads have been joined to
    // ensure the computation has been completed.
    EXPECT_EQ( 31, first_odd );
    EXPECT_EQ( 10, first_even );
}

//---------------------------------------------------------------------------//
//                 end of tstThread.cc
//---------------------------------------------------------------------------//
