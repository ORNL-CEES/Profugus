//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstMatrixVectorMultiply.cc
 * \author Stuart Slattery
 * \brief  HPX matrix-vector multiply test. Effectively an OpenMP style
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
#include <cmath>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Vector.
//---------------------------------------------------------------------------//
class Vector
{
  public:
    using data_it = std::vector<double>::iterator;
    using const_data_it = std::vector<double>::const_iterator;

  public:
    Vector() : d_length( 0 ) { /* ... */ }
    Vector( const int length )
	: d_length( length )
	, d_data( length, 0.0 )
    { /* ... */ }

    int length() const { return d_length; }
    void set( const int i, const double value ) { d_data[i] = value; }
    double get( const int i ) { return d_data[i]; }
    std::vector<double>& raw_data() { return d_data; }
    const std::vector<double>& raw_data() const { return d_data; }

  private:
    int d_length;
    std::vector<double> d_data;
};

//---------------------------------------------------------------------------//
// Matrix
//---------------------------------------------------------------------------//
class Matrix
{
  public:

    Matrix() 
	: d_num_rows( 0 )
	, d_num_cols( 0 )
    { /* ... */ }

    Matrix( const int num_rows, const int num_cols )
	: d_num_rows( num_rows )
	, d_num_cols( num_cols )
	, d_data( num_rows, std::vector<double>(num_cols, 0.0) )
    { /* ... */ }

    int rows() const { return d_num_rows; }
    int cols() const { return d_num_cols; }

    void set( const int row, const int col, const double value )
    { d_data[row][col] = value; }

    double get( const int row, const int col ) const
    { return d_data[row][col]; }

    std::vector<double>& raw_data( const int row )
    { return d_data[row]; }

    const std::vector<double>& raw_data( const int row ) const
    { return d_data[row]; }

  private:

    // Number of rows.
    int d_num_rows;

    // Number of columns.
    int d_num_cols;

    // Local matrix data. Row-major ordering.
    std::vector<std::vector<double> > d_data;
};

//---------------------------------------------------------------------------//
// Opaque task holder. Provides a mechanism to wait on an arbitrary task
// with a future that has a return value we dont care about.
//---------------------------------------------------------------------------//
class OpaqueTaskImplBase
{
  public:
    virtual ~OpaqueTaskImplBase() = default;
    virtual void wait() = 0;
};

//---------------------------------------------------------------------------//
template<class T>
class OpaqueTaskImpl : public OpaqueTaskImplBase
{
  public:
    using future_value_type = T;

  public:    
    OpaqueTaskImpl( hpx::future<T>&& future )
	: d_future( std::move(future) )
    { /* ... */ }

    OpaqueTaskImpl( const OpaqueTaskImpl& opaque ) = delete;
    OpaqueTaskImpl& operator=( const OpaqueTaskImpl& opaque ) = delete;

    OpaqueTaskImpl( OpaqueTaskImpl&& opaque ) = default;
    OpaqueTaskImpl& operator=( OpaqueTaskImpl&& opaque ) = default;

    void wait() override { d_future.wait(); }

  private:
    hpx::future<T> d_future;
};

//---------------------------------------------------------------------------//
class OpaqueTask
{
  public:

    template<class T>
    OpaqueTask( hpx::future<T>&& future )
    {
	d_task_impl = std::unique_ptr<OpaqueTaskImplBase>( 
	    new OpaqueTaskImpl<T>(std::move(future)) );
    }

    OpaqueTask( const OpaqueTask& opaque ) = delete;
    OpaqueTask& operator=( const OpaqueTask& opaque ) = delete;

    OpaqueTask( OpaqueTask&& opaque ) = default;
    OpaqueTask& operator=( OpaqueTask&& opaque ) = default;

    void wait() { d_task_impl->wait(); }

  private:
    std::unique_ptr<OpaqueTaskImplBase> d_task_impl;
};

//---------------------------------------------------------------------------//
// Vector operations.
//---------------------------------------------------------------------------//
class VectorOps
{
  public:

    static void fill_sync( Vector& x, const double a )
    {
	auto fill_op = [&x,a]( const int i ){ x.set(i,a); };
	auto range = boost::irange(0,x.length());
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 fill_op );
    }

    static OpaqueTask fill_async( Vector& x, const double a )
    {
	auto fill_op = [&x,a]( const int i ){ x.set(i,a); };
	auto range = boost::irange(0,x.length());
	return OpaqueTask(
	    hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				     std::begin(range),
				     std::end(range),
				     fill_op )
	    );
    }

    static void scale_sync( Vector& x, const double a )
    {
	auto scale_op = [&x,a]( const int i ){ x.raw_data()[i] *= a; };
	auto range = boost::irange(0,x.length());
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 scale_op );
    }

    static OpaqueTask scale_async( Vector& x, const double a )
    {
	auto scale_op = [&x,a]( const int i ){ x.raw_data()[i] *= a; };
	auto range = boost::irange(0,x.length());
	return OpaqueTask(
	    hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				     std::begin(range),
				     std::end(range),
				     scale_op )
	    );
    }

    static void add_sync( const Vector& x, const Vector& y, Vector& z )
    {
	HPX_ASSERT( z.length() == x.length() );
	HPX_ASSERT( z.length() == y.length() );

	auto sum_op = [&x,&y,&z]( const int i )
		      { z.raw_data()[i] = x.raw_data()[i] + y.raw_data()[i]; };

	auto range = boost::irange(0,x.length());
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 sum_op );
    }

    static OpaqueTask add_async( const Vector& x, const Vector& y, Vector& z )
    {
	HPX_ASSERT( z.length() == x.length() );
	HPX_ASSERT( z.length() == y.length() );

	auto sum_op = [&x,&y,&z]( const int i )
		      { z.raw_data()[i] = x.raw_data()[i] + y.raw_data()[i]; };

	auto range = boost::irange(0,x.length());
	return OpaqueTask(
	    hpx::parallel::for_each( hpx::parallel::parallel_task_execution_policy(),
				     std::begin(range),
				     std::end(range),
				     sum_op )
	    );
    }

    static double dot_sync( const Vector& x, const Vector& y )
    {
	HPX_ASSERT( x.length() == y.length() );
	return hpx::parallel::inner_product( 
	    hpx::parallel::parallel_execution_policy(),
	    std::begin(x.raw_data()),
	    std::end(x.raw_data()),
	    std::begin(y.raw_data()),
	    0.0 );
    }

    static hpx::future<double> dot_async( const Vector& x, const Vector& y )
    {
	HPX_ASSERT( x.length() == y.length() );
	return hpx::parallel::inner_product( 
	    hpx::parallel::parallel_task_execution_policy(),
	    std::begin(x.raw_data()),
	    std::end(x.raw_data()),
	    std::begin(y.raw_data()),
	    0.0 );
    }

    static double square_root( const double a )
    {
	return std::sqrt(a);
    }

    static double norm2_sync( const Vector& x )
    {
	return square_root( dot_sync(x,x) );
    }

    static hpx::future<double> norm2_async( const Vector& x )
    {
	return hpx::lcos::local::dataflow( hpx::launch::async,
					   hpx::util::unwrapped(square_root), 
					   dot_async(x,x) );
    }
};

//---------------------------------------------------------------------------//
// Matrix operations.
//---------------------------------------------------------------------------//
class MatrixOps
{
  public:

    static void apply_sync( const Matrix& A, const Vector& x, Vector& y )
    {
	HPX_ASSERT( A.cols() == x.length()  );
	HPX_ASSERT( A.rows() == y.length()  );

	// Get the vector data.
	const std::vector<double>& x_data = x.raw_data();
	std::vector<double>& y_data = y.raw_data();

	// Create a sum operator for the threads.
	auto sum_op = [&]( const int row )
		      {
			  const std::vector<double>& A_data = A.raw_data( row );
			  y_data[row] = 0.0;
			  for ( auto mij = std::begin(A_data),
				      xj = std::begin(x_data);
				mij != std::end(A_data);
				++mij, ++xj )	 
			      y_data[row] += (*mij) * (*xj);
		      };

	// Loop in parallel over all the rows.
	auto range = boost::irange( 0, A.rows() );
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 sum_op );
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( vector, vector_sync_test )
{
    int N = 100000;
    Vector x( N );
    Vector y( N );

    EXPECT_EQ( N, x.length() );
    EXPECT_EQ( N, y.length() );

    for ( int i = 0; i < N; ++i )
    {
	x.set( i, 1.0 );
	y.set( i, 1.0 );
    }

    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 1.0 );
	EXPECT_EQ( y.get(i), 1.0 );
    }

    EXPECT_EQ( VectorOps::dot_sync(x,y), N );
    EXPECT_EQ( VectorOps::norm2_sync(x), std::sqrt(N) );
    EXPECT_EQ( VectorOps::norm2_sync(y), std::sqrt(N) );

    Vector z( N );
    EXPECT_EQ( N, z.length() );

    VectorOps::add_sync( x, y, z );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( z.get(i), 2.0 );
    }
    EXPECT_EQ( VectorOps::norm2_sync(z), std::sqrt(4*N) );

    VectorOps::fill_sync( x, 3.3 );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 3.3 );
    }

    VectorOps::scale_sync( x, 2.0 );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 6.6 );
    }
}

//---------------------------------------------------------------------------//
TEST( vector, vector_async_test )
{
    int N = 100000;
    Vector x( N );
    Vector y( N );

    EXPECT_EQ( N, x.length() );
    EXPECT_EQ( N, y.length() );

    for ( int i = 0; i < N; ++i )
    {
	x.set( i, 1.0 );
	y.set( i, 1.0 );
    }

    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 1.0 );
	EXPECT_EQ( y.get(i), 1.0 );
    }

    hpx::future<double> xy_dot = VectorOps::dot_async(x,y);
    hpx::future<double> x_norm2 = VectorOps::norm2_async(x);
    hpx::future<double> y_norm2 = VectorOps::norm2_async(y);
    EXPECT_EQ( xy_dot.get(), N );
    EXPECT_EQ( x_norm2.get(), std::sqrt(N) );
    EXPECT_EQ( y_norm2.get(), std::sqrt(N) );

    Vector z( N );
    EXPECT_EQ( N, z.length() );

    OpaqueTask xy_add = VectorOps::add_async( x, y, z );
    xy_add.wait();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( z.get(i), 2.0 );
    }
    hpx::future<double> z_norm2 = VectorOps::norm2_async(z);
    EXPECT_EQ( z_norm2.get(), std::sqrt(4*N) );

    OpaqueTask x_fill = VectorOps::fill_async( x, 3.3 );
    x_fill.wait();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 3.3 );
    }

    OpaqueTask x_scale = VectorOps::scale_async( x, 2.0 );
    x_scale.wait();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 6.6 );
    }
}

//---------------------------------------------------------------------------//
TEST( matrix, matrix_sync_test )
{
    int N = 10000;
    Matrix A( N, N );
    Vector x( N );
    Vector y( N );

    EXPECT_EQ( N, A.rows() );
    EXPECT_EQ( N, A.cols() );

    for ( int i = 0; i < N; ++i )
    {
	A.set( i, i, 1.0 );
	x.set( i, i );
    }

    MatrixOps::apply_sync( A, x, y );
    for ( int i = 0; i < N; ++i )
    {
    	EXPECT_EQ( x.get(i), y.get(i) );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstMatrixVectorMultiply.cc
//---------------------------------------------------------------------------//
