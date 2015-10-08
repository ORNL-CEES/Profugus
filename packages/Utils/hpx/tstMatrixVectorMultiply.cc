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
class OpaqueTaskFutureImplBase
{
  public:
    virtual ~OpaqueTaskFutureImplBase() = default;
    virtual hpx::future<void> get_future() = 0;
};

//---------------------------------------------------------------------------//
template<class T>
class OpaqueTaskFutureImpl : public OpaqueTaskFutureImplBase
{
  public:
    using future_value_type = T;

  public:    
    OpaqueTaskFutureImpl( hpx::future<T>&& future )
	: d_future( std::move(future) )
    { /* ... */ }

    OpaqueTaskFutureImpl( const OpaqueTaskFutureImpl& opaque ) = delete;
    OpaqueTaskFutureImpl& operator=( const OpaqueTaskFutureImpl& opaque ) = delete;

    OpaqueTaskFutureImpl( OpaqueTaskFutureImpl&& opaque ) = default;
    OpaqueTaskFutureImpl& operator=( OpaqueTaskFutureImpl&& opaque ) = default;

    hpx::future<void> get_future() override 
    { 
	return hpx::lcos::local::dataflow( 
	    hpx::launch::async, 
	    hpx::util::unwrapped(OpaqueTaskFutureImpl<T>::flow_trigger), 
	    d_future );
    }

    static void flow_trigger( const T trigger_val ) { /* ... */ }

  private:
    hpx::future<T> d_future;
};

//---------------------------------------------------------------------------//
class OpaqueTaskFuture
{
  public:

    template<class T>
    OpaqueTaskFuture( hpx::future<T>&& future )
    {
	d_task_impl = std::unique_ptr<OpaqueTaskFutureImplBase>( 
	    new OpaqueTaskFutureImpl<T>(std::move(future)) );
    }

    OpaqueTaskFuture( const OpaqueTaskFuture& opaque ) = delete;
    OpaqueTaskFuture& operator=( const OpaqueTaskFuture& opaque ) = delete;

    OpaqueTaskFuture( OpaqueTaskFuture&& opaque ) = default;
    OpaqueTaskFuture& operator=( OpaqueTaskFuture&& opaque ) = default;

    hpx::future<void> get_future() { return d_task_impl->get_future(); }

  private:
    std::unique_ptr<OpaqueTaskFutureImplBase> d_task_impl;
};

//---------------------------------------------------------------------------//
// Operation Tags
//---------------------------------------------------------------------------//
class SyncTag
{
  public:
    template<class T> using value_return_type = T;
    using task_return_type = void;
    using execution_policy = hpx::parallel::parallel_execution_policy;

    template<class T>
    static value_return_type<T> value_return_wrapper( hpx::future<T>&& input ) 
    { return input.get(); }
    
    template<class Input>
    static task_return_type task_return_wrapper( Input&& input ) { /* ... */ }
};

class AsyncTag 
{
  public:
    template<class T> using value_return_type = hpx::future<T>;
    using task_return_type = OpaqueTaskFuture;
    using execution_policy = hpx::parallel::parallel_task_execution_policy;

    template<class T>
    static value_return_type<T>&& value_return_wrapper( hpx::future<T>&& input ) 
    { return std::move(input); }

    template<class Input>
    static task_return_type task_return_wrapper( Input&& input ) 
    { return task_return_type(std::move(input)); }
};

//---------------------------------------------------------------------------//
// Vector operations.
//---------------------------------------------------------------------------//
template<class Tag>
class VectorOps
{
  public:

    // Fill a vector with a scalar.
    static typename Tag::task_return_type fill( Vector& x, const double a ) 
    {
	auto fill_op = [&x,a]( const int i ){ x.set(i,a); };
	auto range = boost::irange(0,x.length());
	return Tag::task_return_wrapper(
	    hpx::parallel::for_each( typename Tag::execution_policy(),
				     std::begin(range),
				     std::end(range),
				     fill_op )
	    );
    }

    // Scale a vector by a scalar.
    static typename Tag::task_return_type scale( Vector& x, const double a )
    {
	auto scale_op = [&x,a]( const int i ){ x.raw_data()[i] *= a; };
	auto range = boost::irange(0,x.length());
	return Tag::task_return_wrapper(
	    hpx::parallel::for_each( typename Tag::execution_policy(),
				     std::begin(range),
				     std::end(range),
				     scale_op )
	    );
    }

    // Add two vectors together and put the results in a third vector.
    static typename Tag::task_return_type 
    add( const Vector& x, const Vector& y, Vector& z )
    {
	HPX_ASSERT( z.length() == x.length() );
	HPX_ASSERT( z.length() == y.length() );
	auto sum_op = [&x,&y,&z]( const int i )
		      { z.raw_data()[i] = x.raw_data()[i] + y.raw_data()[i]; };
	auto range = boost::irange(0,x.length());
	return Tag::task_return_wrapper(
	    hpx::parallel::for_each( typename Tag::execution_policy(),
				     std::begin(range),
				     std::end(range),
				     sum_op )
	    );
    }

    // Compute the dot product of two vectors.
    static typename Tag::template value_return_type<double> 
    dot( const Vector& x, const Vector& y )
    {
	HPX_ASSERT( x.length() == y.length() );
	return hpx::parallel::inner_product( 
	    typename Tag::execution_policy(),
	    std::begin(x.raw_data()),
	    std::end(x.raw_data()),
	    std::begin(y.raw_data()),
	    0.0 );
    }

    // Compute the square root of a number.
    static double square_root( const double a )
    {
	return std::sqrt(a);
    }

    // Compute the 2 norm of a vector.
    static typename Tag::template value_return_type<double> 
    norm2( const Vector& x )
    {
	return Tag::value_return_wrapper(
	    hpx::lcos::local::dataflow( hpx::launch::async,
					hpx::util::unwrapped(square_root), 
					VectorOps<AsyncTag>::dot(x,x) )
	    );
    }
};

//---------------------------------------------------------------------------//
// Matrix operations.
//---------------------------------------------------------------------------//
template<class Tag>
class MatrixOps
{
  public:

    static typename Tag::task_return_type
    apply( const Matrix& A, const Vector& x, Vector& y )
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
	return Tag::task_return_wrapper(
	    hpx::parallel::for_each( typename Tag::execution_policy(),
				     std::begin(range),
				     std::end(range),
				     sum_op )
	    );
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

    EXPECT_EQ( VectorOps<SyncTag>::dot(x,y), N );
    EXPECT_EQ( VectorOps<SyncTag>::norm2(x), std::sqrt(N) );
    EXPECT_EQ( VectorOps<SyncTag>::norm2(y), std::sqrt(N) );

    Vector z( N );
    EXPECT_EQ( N, z.length() );

    VectorOps<SyncTag>::add( x, y, z );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( z.get(i), 2.0 );
    }
    EXPECT_EQ( VectorOps<SyncTag>::norm2(z), std::sqrt(4*N) );

    VectorOps<SyncTag>::fill( x, 3.3 );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 3.3 );
    }

    VectorOps<SyncTag>::scale( x, 2.0 );
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

    MatrixOps<SyncTag>::apply( A, x, y );
    for ( int i = 0; i < N; ++i )
    {
    	EXPECT_EQ( x.get(i), y.get(i) );
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

    Vector z( N );
    EXPECT_EQ( N, z.length() );

    auto xy_dot = VectorOps<AsyncTag>::dot(x,y);
    auto x_norm2 = VectorOps<AsyncTag>::norm2(x);
    auto y_norm2 = VectorOps<AsyncTag>::norm2(y);
    auto xy_add = VectorOps<AsyncTag>::add( x, y, z );

    EXPECT_EQ( xy_dot.get(), N );
    EXPECT_EQ( x_norm2.get(), std::sqrt(N) );
    EXPECT_EQ( y_norm2.get(), std::sqrt(N) );

    xy_add.get_future().wait();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( z.get(i), 2.0 );
    }
    auto z_norm2 = VectorOps<AsyncTag>::norm2(z);
    EXPECT_EQ( z_norm2.get(), std::sqrt(4*N) );

    auto x_fill = VectorOps<AsyncTag>::fill( x, 3.3 );
    x_fill.get_future().wait();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 3.3 );
    }

    auto x_scale = VectorOps<AsyncTag>::scale( x, 2.0 );
    x_scale.get_future().wait();
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 6.6 );
    }
}

//---------------------------------------------------------------------------//
TEST( matrix, matrix_async_test )
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

    auto apply_op = MatrixOps<AsyncTag>::apply( A, x, y );
    apply_op.get_future().wait();
    for ( int i = 0; i < N; ++i )
    {
    	EXPECT_EQ( x.get(i), y.get(i) );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstMatrixVectorMultiply.cc
//---------------------------------------------------------------------------//
