//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstConditionVariable.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 9:49:03 2015
 * \brief  condition variable testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// HELPER DATA
//---------------------------------------------------------------------------//
class TestData
{
  public:

    // Constructor.
    TestData()
	: d_ready( false )
    { /* ... */ }

    // Set the data.
    void set_data( std::vector<int>&& data )
    {
	d_x = data;
    }

    // Find the first even value in the class data.
    void even_func( int& i )
    {
	wait_and_unlock();
	i = *std::find_if(std::begin(d_x), std::end(d_x),
			  [](int n){ return n % 2 == 0; });
    }

    // Find the first odd value in the class data.
    void odd_func( int& i )
    {
	wait_and_unlock();
	i = *std::find_if(std::begin(d_x), std::end(d_x),
			  [](int n){ return n % 2 == 1; });
    }

    // Execute all functions that are waiting.
    void run_funcs()
    {
	std::lock_guard<std::mutex> guard( d_mutex );
	d_ready = true;
	d_run_condition.notify_all();
    }

  private:

    // Wait for the ready signal.
    void wait_and_unlock()
    {
	std::unique_lock<std::mutex> lock( d_mutex );
	d_run_condition.wait( lock, [this](){return d_ready;} );
	lock.unlock();
    }

  private:

    // Class data.
    std::vector<int> d_x;

    // Data mutex.
    std::mutex d_mutex;

    // Execution condition boolean.
    bool d_ready;

    // Execution condition variable.
    std::condition_variable d_run_condition;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(condition_variable, condition_variable_test)
{
    // Create some test data.
    TestData test_data;
    test_data.set_data( {10, 22, 31, 44, 56} );

    // initialize a thread to find first odd
    int first_odd = -1;
    std::thread odd_thread( 
	&TestData::odd_func, &test_data, std::ref(first_odd) );

    // initialize a thread to find first even
    int first_even = -1;
    std::thread even_thread( 
	&TestData::even_func, &test_data, std::ref(first_even) );
    
    // check that the threads have different ids
    EXPECT_NE( odd_thread.get_id(), even_thread.get_id() );

    // run the threads
    test_data.run_funcs();
    
    // join the threads before before checking the test results
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
//                 end of tstConditionVariable.cc
//---------------------------------------------------------------------------//
