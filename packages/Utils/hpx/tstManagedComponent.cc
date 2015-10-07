//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstManagedComponent.cc
 * \author Stuart Slattery
 * \brief  HPX managed component testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <hpx/hpx.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/components.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>
#include <hpx/runtime/components/component_factory.hpp>
#include <hpx/runtime/components/stubs/stub_base.hpp>
#include <hpx/runtime/applier/apply.hpp>

#include <mutex>
#include <vector>
#include <algorithm>
#include <string>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
/* 
   Managed components are basically data structures where hpx manages all of
   the locking automatically as well as the location of the data once an
   inheritance structure is implemented. The data structure is then thread
   safe with locks. This component is managed because we might want to create
   it in bulk (i.e. a lot of them). A server where the data actually resides
   and can be on any locality while a client provides an interface to the data
   on a locality (but could be called from elsewhere).
 */

//---------------------------------------------------------------------------//
// Test Helpers
//---------------------------------------------------------------------------//
// Server class. This is the basic class that actually does the work. If
// refactoring existing code to a make thread-safe, this would be the starting
// point. We also add a 'testing' namespace here. Otherwise the server
// namespace will clash with hpx namespaces in the macros.
namespace testing
{
namespace server
{
class DataBank : public hpx::components::locking_hook<
    hpx::components::managed_component_base<DataBank> >
{
  public:
    DataBank() : d_data( 0.0 ) { /* ... */ }
    DataBank( const double init_value ) : d_data( init_value ) { /* ... */ }

    void reset() { d_data = 0.0; }
    void add( const double value ) { d_data += value; }
    double get() { return d_data; }

    HPX_DEFINE_COMPONENT_ACTION( DataBank, reset );
    HPX_DEFINE_COMPONENT_ACTION( DataBank, add );
    HPX_DEFINE_COMPONENT_ACTION( DataBank, get );

  private:
    double d_data;
};
} // end namespace server
} // end namespace testing

// Now declare the action registration the were created from the component
// definition. HPX will add an _action suffix to those functions it defined in
// the above macro.
HPX_REGISTER_ACTION_DECLARATION( testing::server::DataBank::reset_action,
				 data_bank_reset_action );
HPX_REGISTER_ACTION_DECLARATION( testing::server::DataBank::add_action,
				 data_bank_add_action );
HPX_REGISTER_ACTION_DECLARATION( testing::server::DataBank::get_action,
				 data_bank_get_action );

// Now register the component with the factory.
HPX_REGISTER_COMPONENT_MODULE();
HPX_REGISTER_COMPONENT( 
    hpx::components::managed_component<testing::server::DataBank>,
    data_bank );

// Finally, register the actions.
HPX_REGISTER_ACTION( testing::server::DataBank::reset_action,
		     data_bank_reset_action );
HPX_REGISTER_ACTION( testing::server::DataBank::add_action,
		     data_bank_add_action );
HPX_REGISTER_ACTION( testing::server::DataBank::get_action,
		     data_bank_get_action );

//---------------------------------------------------------------------------//
// Stubs class. Stateless class that manages the execution of actions
// (functions) on the server class for a given target. These are helper
// functions and could technically be bypassed and implemented directly in the
// client.
namespace testing
{
namespace stubs
{
class DataBank : public hpx::components::stub_base<server::DataBank>
{
  public:
    // Non-blocking version of reset. Using hpx::apply instead of async
    // produces a non-blocking version (basically a thread launch without a
    // future and then detaching the thread). We can do this because reset()
    // returns void.
    static void reset_non_blocking( const hpx::naming::id_type& gid )
    {
	hpx::apply<server::DataBank::reset_action>( gid );
    }

    // Blocking version where one has to wait on a returned future before
    // exiting. 
    static void reset_sync( const hpx::naming::id_type& gid )
    {
	hpx::async<server::DataBank::reset_action>( gid ).get();
    }

    // Non-blocking add.
    static void add_non_blocking( const hpx::naming::id_type& gid,
				  const double value )
    {
	hpx::apply<server::DataBank::add_action>( gid, value );
    }

    // Blocking add.
    static void add_sync( const hpx::naming::id_type& gid,
			  const double value )
    {
	hpx::async<server::DataBank::add_action>( gid, value ).get();
    }
    
    // Asynchronous get function. We can't do a non-blocking because this
    // function returns a value. Instead we return a future so we can wait
    // later instead of now.
    static hpx::future<double> get_async( const hpx::naming::id_type& gid )
    {
	return hpx::async<server::DataBank::get_action>( gid );
    }

    // Synchronous get function.
    static double get_sync( const hpx::naming::id_type& gid )
    {
	return hpx::async<server::DataBank::get_action>( gid ).get();
    }
};
} // end namespace stubs
} // end namespace testing

//---------------------------------------------------------------------------//
// Client class. The client is the high-level facade for the thread-safe
// server. It uses the stubs class to manage the server and presents the user
// with a standard API while under the hood the threading is managed.
namespace testing
{
class DataBank : public hpx::components::client_base<DataBank,stubs::DataBank>
{
  public:

    using base_type = hpx::components::client_base<DataBank,stubs::DataBank>;

    // Default constructor. No component is associated with this.
    DataBank()
    { /* ... */ }

    // Create a client for an existing server with a given GID.
    DataBank( hpx::future<hpx::naming::id_type>&& gid )
	: base_type( std::move(gid) )
    { /* ... */ }

    // Client interface for the DataBank server.
    void reset_non_blocking()
    {
	HPX_ASSERT( this->get_id() );
	this->base_type::reset_non_blocking( this->get_id() );
    }

    void reset_sync()
    {
	HPX_ASSERT( this->get_id() );
	this->base_type::reset_sync( this->get_id() );
    }

    void add_non_blocking( const double value )
    {
	HPX_ASSERT( this->get_id() );
	this->base_type::add_non_blocking( this->get_id(), value );
    }

    void add_sync( const double value)
    {
	HPX_ASSERT( this->get_id() );
	this->base_type::add_sync( this->get_id(), value );
    }

    hpx::future<double> get_async()
    {
	HPX_ASSERT( this->get_id() );
	return this->base_type::get_async( this->get_id() );
    }

    double get_sync()
    {
	HPX_ASSERT( this->get_id() );
	return this->base_type::get_sync( this->get_id() );
    }
};
} // end namespace testing

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( managed_component, add_numbers_test )
{
    // Create a data bank client with a local server.
    testing::DataBank bank( 
	hpx::components::new_<testing::server::DataBank>(hpx::find_here())
	);

    // Test the bank.
    EXPECT_EQ( bank.get_sync(), 0.0 );

    bank.add_non_blocking( 3.3 );
    EXPECT_EQ( bank.get_sync(), 3.3 );    

    hpx::future<double> bf = bank.get_async();
    bank.add_sync( 4.0 );
    EXPECT_EQ( bf.get(), 3.3 );
    EXPECT_EQ( bank.get_sync(), 7.3 );

    bank.reset_non_blocking();
    EXPECT_EQ( bank.get_sync(), 0.0 );

    bank.add_non_blocking( 2.2 );
    EXPECT_EQ( bank.get_sync(), 2.2 );

    bf = bank.get_async();
    bank.reset_sync();
    EXPECT_EQ( bf.get(), 2.2 );
    EXPECT_EQ( bank.get_sync(), 0.0 );
}

//---------------------------------------------------------------------------//
//                 end of tstManagedComponent.cc
//---------------------------------------------------------------------------//
