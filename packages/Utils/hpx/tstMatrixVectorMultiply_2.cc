//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstMatrixVectorMultiply.cc
 * \author Stuart Slattery
 * \brief  HPX matrix-vector multiply test.
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
// Vector operations.
//---------------------------------------------------------------------------//
class VectorOps
{
  public:

    static void fill( Vector& x, const double a )
    {
	auto fill_op = [&x=x,&a]( const int i ){ x.set(i,a); };
	auto range = boost::irange(0,x.length());
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 fill_op );
    }

    static void scale( Vector& x, const double a )
    {
	auto scale_op = [&x=x,&a]( const int i ){ x.raw_data()[i] *= a; };
	auto range = boost::irange(0,x.length());
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 scale_op );
    }

    static void add( const Vector& x, const Vector& y, Vector& z )
    {
	HPX_ASSERT( z.length() == x.length() );
	HPX_ASSERT( z.length() == y.length() );

	auto sum_op = [&x,&y,&z=z]( const int i )
		      { z.raw_data()[i] = x.raw_data()[i] + y.raw_data()[i]; };

	auto range = boost::irange(0,x.length());
	hpx::parallel::for_each( hpx::parallel::parallel_execution_policy(),
				 std::begin(range),
				 std::end(range),
				 sum_op );
    }

    static double dot( const Vector& x, const Vector& y )
    {
	HPX_ASSERT( x.length() == y.length() );
	return hpx::parallel::inner_product( 
	    hpx::parallel::parallel_execution_policy(),
	    std::begin(x.raw_data()),
	    std::end(x.raw_data()),
	    std::begin(y.raw_data()),
	    0.0 );
    }

    static double norm2( const Vector& x )
    {
	return std::sqrt( dot(x,x) );
    }
};

//---------------------------------------------------------------------------//
// Matrix operations.
//---------------------------------------------------------------------------//
class MatrixOps
{
  public:

    static void apply( const Matrix& A, const Vector& x, Vector& y )
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
TEST( vector, vector_test )
{
    int N = 10;
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

    EXPECT_EQ( VectorOps::dot(x,y), N );
    EXPECT_EQ( VectorOps::norm2(x), std::sqrt(N) );
    EXPECT_EQ( VectorOps::norm2(y), std::sqrt(N) );

    Vector z( N );
    EXPECT_EQ( N, z.length() );

    VectorOps::add( x, y, z );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( z.get(i), 2.0 );
    }
    EXPECT_EQ( VectorOps::norm2(z), std::sqrt(4*N) );

    VectorOps::fill( x, 3.3 );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 3.3 );
    }

    VectorOps::scale( x, 2.0 );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i), 6.6 );
    }
}

//---------------------------------------------------------------------------//
TEST( matrix, matrix_test )
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

    MatrixOps::apply( A, x, y );
    for ( int i = 0; i < N; ++i )
    {
    	EXPECT_EQ( x.get(i), y.get(i) );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstMatrixVectorMultiply.cc
//---------------------------------------------------------------------------//
