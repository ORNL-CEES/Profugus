//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstAtomic.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 9:49:03 2015
 * \brief  Atomic testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <atomic>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(atomic, atomic_test)
{
    std::atomic<int> x( 0 );

    int last_int = x.fetch_add( 3 );
    EXPECT_EQ( 0, last_int );

    int current_int = x.load();
    EXPECT_EQ( 3, current_int );

    last_int = x.exchange( 23 );
    EXPECT_EQ( 3, last_int );

    current_int = x.load();
    EXPECT_EQ( 23, current_int );

    int expected_int = 23;
    bool did_exchange = x.compare_exchange_strong( expected_int, 44 );
    EXPECT_TRUE( did_exchange );
    EXPECT_EQ( expected_int, 23 );
    
    current_int = x.load();
    EXPECT_EQ( 44, current_int );

    did_exchange = x.compare_exchange_strong( expected_int, 3 );
    EXPECT_FALSE( did_exchange );
    EXPECT_EQ( expected_int, 44 );    
}

//---------------------------------------------------------------------------//
//                 end of tstAtomic.cc
//---------------------------------------------------------------------------//
