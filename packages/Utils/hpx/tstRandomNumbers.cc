//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstRandomNumbers.cc
 * \author Stuart Slattery
 * \brief  HPX conjugate gradient test.
 * implementation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <hpx/hpx.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/components.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>
#include <hpx/runtime/components/component_factory.hpp>
#include <hpx/runtime/components/stubs/stub_base.hpp>
#include <hpx/runtime/applier/apply.hpp>
#include <hpx/parallel/algorithms/inner_product.hpp>
#include <hpx/parallel/algorithms/reduce.hpp>
#include <hpx/parallel/execution_policy.hpp>

#include <boost/range/irange.hpp>

#include <mutex>
#include <vector>
#include <algorithm>
#include <string>
#include <random>
#include <cmath>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// XORSHIFT
//---------------------------------------------------------------------------//
// Server
class Xorshift : public hpx::components::simple_component_base<Xorshift>
{
  public:

    //@{
    //! Typedefs.
    typedef std::uniform_int_distribution<int> uniform_int_distribution_type;
    typedef std::uniform_real_distribution<double> uniform_real_distribution_type;
    typedef uint_fast64_t result_type;
    //@}

    // Default constructor.
    Xorshift() { /* ... */ }

    //! Constructor.
    explicit Xorshift( const result_type seed )
	: d_x( seed )
    { /* ... */ }

    //! Minimum value.
    static result_type min ()
    { return std::numeric_limits<result_type>::min(); }

    //! Maximum value.
    static result_type max ()
    { return std::numeric_limits<result_type>::max(); }

    //! Seed the engine.
    void seed( const result_type seed )
    { d_x = seed; }

    // Get a random number.
    result_type get()
    {
	d_x ^= d_x << 13;
	d_x ^= d_x >> 7;
	d_x ^= d_x << 17;
	return d_x;
    }

    // Register functions.
    HPX_DEFINE_COMPONENT_DIRECT_ACTION( Xorshift, get, xorshift_get_action );

  private:

    // Random number state.
    result_type d_x;
};

// Register component and actions.
HPX_REGISTER_COMPONENT( hpx::components::simple_component<Xorshift>, Xorshift );
HPX_REGISTER_ACTION( Xorshift::xorshift_get_action );

//---------------------------------------------------------------------------//
// Client.
class XorshiftClient : 
    public hpx::components::client_base<XorshiftClient,Xorshift>
{
  public:

    using base_type = hpx::components::client_base<XorshiftClient,Xorshift>;

  public:
    
    XorshiftClient() { /* ... */ }

    XorshiftClient( hpx::id_type where, const Xorshift::result_type seed )
	: base_type( hpx::new_<Xorshift>(where,seed) )
    { /* ... */ }

    XorshiftClient( hpx::id_type where, Xorshift&& rng )
      : base_type( hpx::new_<Xorshift>(hpx::colocated(where),std::move(rng)) )
    { /* ... */ }

    XorshiftClient( hpx::future<hpx::id_type>&& id )
	: base_type( std::move(id) )
    { /* ... */ }

    XorshiftClient( hpx::future<XorshiftClient>&& client )
	: base_type( std::move(client) )
    { /* ... */ }

    hpx::future<Xorshift::result_type> get_async()
    {
	HPX_ASSERT( this->get_id() );
	return hpx::async( Xorshift::xorshift_get_action(), this->get_id() );
    }

    Xorshift::result_type get_sync()
    {
	HPX_ASSERT( this->get_id() );
	return hpx::async( Xorshift::xorshift_get_action(), this->get_id() ).get();
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( xorshift, rng_test )
{
    using Result = Xorshift::result_type;
    Result seed = 39439237;

    XorshiftClient rng_1( hpx::find_here(), seed );
    XorshiftClient rng_2( hpx::find_here(), seed );

    int num_ran = 3;
    std::vector<Result> r1( num_ran );
    std::vector<Result> r2( num_ran );

    r1[0] = rng_1.get_sync();
    r1[1] = rng_1.get_sync();
    r1[2] = rng_1.get_sync();

    auto n0 = rng_2.get_async();
    auto n1 = rng_2.get_async();
    auto n2 = rng_2.get_async();

    r2[2] = n2.get();
    r2[1] = n1.get();
    r2[0] = n0.get();

    EXPECT_EQ( r1[0], r2[0] );
    EXPECT_EQ( r1[1], r2[1] );
    EXPECT_EQ( r1[2], r2[2] );
}

//---------------------------------------------------------------------------//
TEST( xorshift, sync_array_test_1 )
{
    using Result = Xorshift::result_type;
    Result seed = 39439237;

    XorshiftClient rng( hpx::find_here(), seed );

    int num_ran = 1000000;
    std::vector<Result> ran_vec( num_ran );

    auto fill_op = [&](int i){ ran_vec[i] = rng.get_sync(); };

    auto range = boost::irange( 0, num_ran );
    hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
			     std::begin(range),
			     std::end(range),
			     fill_op );
}

//---------------------------------------------------------------------------//
TEST( xorshift, sync_array_test_2 )
{
    using Result = Xorshift::result_type;
    Result seed = 39439237;

    XorshiftClient rng( hpx::find_here(), seed );

    int num_ran = 1000000;
    std::vector<Result> ran_vec( num_ran );

    auto fill_op = [&](int i){ ran_vec[i] = rng.get_sync(); };

    auto range = boost::irange( 0, num_ran );
    auto fill_task =
	hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 fill_op );
    fill_task.wait();    
}

//---------------------------------------------------------------------------//
TEST( xorshift, async_array_test_1 )
{
    using Result = Xorshift::result_type;
    Result seed = 39439237;

    XorshiftClient rng( hpx::find_here(), seed );

    int num_ran = 1000000;
    std::vector<hpx::future<Result> > ran_vec( num_ran );

    auto fill_op = [&](int i){ ran_vec[i] = rng.get_async(); };

    auto range = boost::irange( 0, num_ran );
    hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
			     std::begin(range),
			     std::end(range),
			     fill_op );
    hpx::lcos::wait_all( ran_vec );
}

//---------------------------------------------------------------------------//
TEST( xorshift, async_array_test_2 )
{
    using Result = Xorshift::result_type;
    Result seed = 39439237;

    XorshiftClient rng( hpx::find_here(), seed );

    int num_ran = 1000000;
    std::vector<hpx::future<Result> > ran_vec( num_ran );

    auto fill_op = [&](int i){ ran_vec[i] = rng.get_async(); };

    auto range = boost::irange( 0, num_ran );
    auto fill_task =
	hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 fill_op );
    fill_task.wait();    
    hpx::lcos::wait_all( ran_vec );
}

//---------------------------------------------------------------------------//
//                 end of tstRandomNumbers.cc
//---------------------------------------------------------------------------//
