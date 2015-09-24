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

TEST(Thread, function1)
{
    
    std::vector<int> x = {10, 22, 31, 44, 56};

    // run a thread to find first odd
    auto odd_func = [&]( std::vector<int>& vec )
		    {
			std::find_if(std::begin(vec), std::end(vec),
				     [](int n){ return n % 2 == 1; });
		    };
    std::thread odd_thread( odd_func, x );
    EXPECT_EQ(31, *i);

    // run a thread to find first even
    auto even_func = [&] ( std::vector<int>& vec )
		     {
			 std::find_if(std::begin(vec), std::end(vec),
				      [](int n){ return n % 2 == 0; });
		     }
    std::thread even_thread( even_func, x );
    EXPECT_EQ(10, *i);

    // check that the threads have different ids
    EXPECT_NEQ( odd_thread.get_id(), even_thread.get_id() );
    
    // join the threads before exiting
    EXPECT_TRUE( odd_thread.joinable() );
    EXPECT_TRUE( even_thread.joinable() );
    odd_thread.join();
    even_thread.join();
    EXPECT_FALSE( odd_thread.joinable() );
    EXPECT_FALSE( even_thread.joinable() );
}

//---------------------------------------------------------------------------//
//                 end of tstThread.cc
//---------------------------------------------------------------------------//
