//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstMutex.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 9:49:03 2015
 * \brief  Mutex testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// HELPER DATA
//---------------------------------------------------------------------------//
class TestData
{
  public:

    // Set the data.
    void set_data( std::vector<int>&& data )
    {
	d_x = data;
    }

    // Find the first even value in the class data.
    void even_func( int& i )
    {
	std::lock_guard<std::mutex> lock( d_mutex );
	i = *std::find_if(std::begin(d_x), std::end(d_x),
			  [](int n){ return n % 2 == 0; });
    }

    // Find the first odd value in the class data.
    void odd_func( int& i )
    {
	std::lock_guard<std::mutex> lock( d_mutex );
	i = *std::find_if(std::begin(d_x), std::end(d_x),
			  [](int n){ return n % 2 == 1; });
    }

  private:

    // Class data.
    std::vector<int> d_x;

    // Data mutex.
    std::mutex d_mutex;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(mutex, mutex_test)
{
    // Create some test data.
    TestData test_data;
    test_data.set_data( {10, 22, 31, 44, 56} );

    // run a thread to find first odd
    int first_odd = -1;
    std::thread odd_thread( 
	&TestData::odd_func, &test_data, std::ref(first_odd) );

    // run a thread to find first even
    int first_even = -1;
    std::thread even_thread( 
	&TestData::even_func, &test_data, std::ref(first_even) );
    
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
    EXPECT_EQ(31, first_odd);
    EXPECT_EQ(10, first_even );
}

//---------------------------------------------------------------------------//
//                 end of tstThread.cc
//---------------------------------------------------------------------------//
