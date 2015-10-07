//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstMatrixVectorMultiply.cc
 * \author Stuart Slattery
 * \brief  HPX matrix-vector multiply test.
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
#include <cmath>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Matrix
//---------------------------------------------------------------------------//
class Matrix
{
  public:

    // Default constructor.
    Matrix() 
	: d_num_rows( 0 )
	, d_num_cols( 0 )
    { /* ... */ }

    // Size constructor.
    Matrix( const int num_rows, const int num_cols )
	: d_num_rows( num_rows )
	, d_num_cols( num_cols )
	, d_data( num_cols, std::vector<double>(num_rows, 0.0) )
    { /* ... */ }

    // Get the number of rows.
    int rows() const { return d_num_rows; }

    // Get the number of columns.
    int cols() const { return d_num_cols; }

    // Set a specific value.
    void set( const int row, const int col, const double value )
    {
	d_data[col][row] = value;
    }

    // Get a specific value.
    double get( const int row, const int col ) const
    {
	return d_data[col][row];
    }

    // Add. this = x + y
    void add( const Matrix& x, const Matrix& y )
    {
	HPX_ASSERT( d_num_rows == x.d_num_rows );
	HPX_ASSERT( d_num_cols == x.d_num_cols );
	HPX_ASSERT( d_num_rows == y.d_num_rows );
	HPX_ASSERT( d_num_cols == y.d_num_cols );
	for ( int j = 0; j < d_num_cols; ++j )
	{
	    for ( int i = 0; i < d_num_rows; ++i )
	    {
		d_data[j][i] = x.d_data[j][i] + y.d_data[j][i];
	    }
	}
    }

    // Dot product. Vector-vector overload.
    double dot( const Matrix& x ) const
    {
	HPX_ASSERT( d_num_rows == x.d_num_rows  );
	HPX_ASSERT( d_num_cols == 1 );
	HPX_ASSERT( x.d_num_cols == 1 );
	double sum = 0.0;
	for ( int i = 0; i < d_num_rows; ++i )
	{
	    sum += d_data[0][i] * x.d_data[0][i];
	}
	return sum;
    }

    // Dot product. Matrix-vector overload. y = this*x
    void dot( const Matrix& x, Matrix& y ) const
    {
	HPX_ASSERT( d_num_rows == y.d_num_rows  );
	HPX_ASSERT( d_num_cols == x.d_num_rows  );
	HPX_ASSERT( x.d_num_cols == y.d_num_cols );

	// Create a vector of futures.
	std::vector<hpx::future<void> > sum_futures;
	sum_futures.reserve( d_num_rows );

	// Create a sum operator for the threads.
	auto sum_op = [this,&x,&y=y]( const int row, const int col )
	{
	    y.d_data[col][row] = 0.0;
	    for ( int j = 0; j < this->d_num_cols; ++j ) 
		y.d_data[col][row] += this->d_data[j][row] * x.d_data[col][j];

	};

	// Calculate the sum. Synchronize on each vector.
	for ( int v = 0; v < x.d_num_cols; ++v )
	{
	    // Launch the sums. One task per row.
	    for ( int i = 0; i < d_num_rows; ++i )
	    {
		sum_futures.push_back( hpx::async(sum_op, i, v) );
	    }

	    // Wait on the sums.
	    hpx::lcos::wait_all(sum_futures);
	}
    }

    // 2-norm of a single column matrix (vector).
    double norm2() const
    {
	HPX_ASSERT( d_num_cols == 1 );
	double norm = 0.0;
	for ( auto& d : d_data[0] ) norm += d*d;
	return std::sqrt(norm);
    }

  private:

    // Number of rows.
    int d_num_rows;

    // Number of columns.
    int d_num_cols;

    // Local matrix data. Column-major ordering.
    std::vector<std::vector<double> > d_data;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST( matrix_vector, vector_dot_test )
{
    int N = 10000;
    Matrix x( N, 1 );
    Matrix y( N, 1 );

    for ( int i = 0; i < N; ++i )
    {
	x.set( i, 0, 1.0 );
	y.set( i, 0, 1.0 );
    }

    EXPECT_EQ( x.dot(y), N );
    EXPECT_EQ( x.norm2(), std::sqrt(N) );
}

//---------------------------------------------------------------------------//
TEST( matrix_vector, multiply_test )
{
    int N = 10000;
    Matrix A( N, N );
    Matrix x( N, 1 );
    Matrix y( N, 1 );

    for ( int i = 0; i < N; ++i )
    {
	A.set( i, i, 1.0 );
	x.set( i, 0, i );
    }

    A.dot( x, y );
    for ( int i = 0; i < N; ++i )
    {
	EXPECT_EQ( x.get(i,0), y.get(i,0) );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstMatrixVectorMultiply.cc
//---------------------------------------------------------------------------//
