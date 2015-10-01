//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstFuture.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 9:49:03 2015
 * \brief  Future testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <exception>
#include <vector>
#include <algorithm>
#include <functional>
#include <future>

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
    int even_func()
    {
	return *std::find_if(std::begin(d_x), std::end(d_x),
			     [](int n){ return n % 2 == 0; });
    }

    // Find the first odd value in the class data.
    int odd_func()
    {
	return *std::find_if(std::begin(d_x), std::end(d_x),
			     [](int n){ return n % 2 == 1; });
    }

  private:

    // Class data.
    std::vector<int> d_x;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(future, future_test)
{
    // Create some test data.
    TestData test_data;
    test_data.set_data( {10, 22, 31, 44, 56} );

    // create a future for the first odd.
    std::future<int> first_odd = std::async( &TestData::odd_func, &test_data );

    // create a future for the first even.
    std::future<int> first_even = std::async( &TestData::even_func, &test_data );

    // check the even and odd results
    EXPECT_EQ(31, first_odd.get() );
    EXPECT_EQ(10, first_even.get() );
}

//---------------------------------------------------------------------------//
TEST( packaged_task, packaged_task_test )
{
    // Create some test data.
    TestData test_data;
    test_data.set_data( {10, 22, 31, 44, 56} );

    // create a packaged task for the first odd.
    std::packaged_task<int()> first_odd( 
	std::bind(&TestData::odd_func, &test_data) );

    // create a packaged task for the first even.
    std::packaged_task<int()> first_even( 
	std::bind(&TestData::even_func, &test_data) );

    // get the task futures
    std::future<int> odd_future = first_odd.get_future();
    std::future<int> even_future = first_even.get_future();

    // execute the tasks by moving the tasks into threads
    std::thread odd_thread( std::move(first_odd) );
    std::thread even_thread( std::move(first_even) );

    // detach the threads. we know they will finish because we will wait on
    // the futures
    odd_thread.detach();
    even_thread.detach();

    // check the even and odd results
    EXPECT_EQ(31, odd_future.get() );
    EXPECT_EQ(10, even_future.get() );    
}

//---------------------------------------------------------------------------//
TEST( promise, promise_test )
{
    // Create a data promise
    std::promise<std::vector<int> > data_promise;

    // Create a shared future for the data promise so we can share the data
    // with multiple tasks
    std::shared_future<std::vector<int> > 
	data_future( data_promise.get_future() );

    // create a task to find the first odd
    auto odd_func = []( std::shared_future<std::vector<int> >& data_future )
		    {
			return *std::find_if(std::begin(data_future.get()), 
					     std::end(data_future.get()),
					     [](int n){ return n % 2 == 1; });
		    };
    std::future<int> first_odd = 
	std::async( odd_func, std::ref(data_future) );

    // create a task to find the first even
    auto even_func = []( std::shared_future<std::vector<int> >& data_future )
		     {
			 return *std::find_if(std::begin(data_future.get()), 
					      std::end(data_future.get()),
					      [](int n){ return n % 2 == 0; });
		     };
    std::future<int> first_even = 
	std::async( even_func, std::ref(data_future) );

    // set the data with the promise
    data_promise.set_value( {10, 22, 31, 44, 56} );

    // check the even and odd results
    EXPECT_EQ(31, first_odd.get() );
    EXPECT_EQ(10, first_even.get() );
}

//---------------------------------------------------------------------------//
TEST( future, future_exception_test )
{
    // this function throws if anything but zero is provided
    auto i_want_zero = []( const int number )
		       {
			   if ( number !=  0 )
			   {
			       throw std::logic_error("I wanted zero!");
			   }
			   return number;
		       };

    // launch a good task
    std::future<int> good_future = std::async( i_want_zero, 0 );

    // launch a bad task
    std::future<int> bad_future = std::async( i_want_zero, 42 );

    // check the good future
    EXPECT_EQ( 0, good_future.get() );

    // make sure we catch a std::logic_error exception when we check the bad
    // future
    bool exception_was_caught = false;
    std::string exception_message;
    try
    {
	int bad_result = bad_future.get();
	EXPECT_EQ( 42, bad_result );
	exception_was_caught = false;
    }
    catch( const std::logic_error& e )
    {
	exception_was_caught = true;
	exception_message = e.what();
    }
    EXPECT_TRUE( exception_was_caught );    
    EXPECT_EQ( "I wanted zero!", exception_message );
}

//---------------------------------------------------------------------------//
//                 end of tstFuture.cc
//---------------------------------------------------------------------------//
