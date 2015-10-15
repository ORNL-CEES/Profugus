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

#include <atomic>
#include <mutex>
#include <vector>
#include <algorithm>
#include <string>
#include <random>
#include <cmath>
#include <map>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
namespace Testing
{
class Xorshift
{
  public:

    //@{
    //! Typedefs.
    typedef std::uniform_int_distribution<int> uniform_int_distribution_type;
    typedef std::uniform_real_distribution<double> uniform_real_distribution_type;
    typedef uint_fast64_t result_type;
    //@}

    //! Minimum value.
    static result_type min()
    { return std::numeric_limits<result_type>::min(); }

    //! Maximum value.
    static result_type max()
    { return std::numeric_limits<result_type>::max(); }

    // Get a random number.
    static void update_state( result_type& x )
    {
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
    }
};
} // end namespace Testing.

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( xorshift, array_test_1 )
{
    using Result = Testing::Xorshift::result_type;
    Result seed = 39439237;

    int num_ran = 10000000;
    std::vector<Result> ran_vec( num_ran, seed );

    auto fill_op = [&](int i){ Testing::Xorshift::update_state(ran_vec[i]); };

    auto range = boost::irange( 0, num_ran );
    hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
			     std::begin(range),
			     std::end(range),
			     fill_op );
}

//---------------------------------------------------------------------------//
TEST( xorshift, array_test_2 )
{
    using Result = Testing::Xorshift::result_type;
    Result seed = 39439237;

    int num_ran = 10000000;
    std::vector<Result> ran_vec( num_ran, seed );

    auto fill_op = [&](int i){ Testing::Xorshift::update_state(ran_vec[i]); };

    auto range = boost::irange( 0, num_ran );
    auto fill_task =
	hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 fill_op );
    fill_task.wait();    
}

//---------------------------------------------------------------------------//
//                 end of tstRandomNumbers.cc
//---------------------------------------------------------------------------//
