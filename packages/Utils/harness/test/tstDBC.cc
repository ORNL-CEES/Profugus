//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstDBC.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 03 20:30:12 2013
 * \brief  DBC unit teset.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "../DBC.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Catch, runtime_error)
{
    bool caught = false;
    try
    {
        throw std::runtime_error( "hello1" );
    }
    catch( const profugus::assertion &a )
    {
        caught = false;
    }
    catch( ... )
    {
        caught = true;
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//

TEST(Catch, profugus_assertion)
{
    bool caught = false;

    std::string error_message;
    try
    {
        throw profugus::assertion( "hello1", "myfile", 42 );
    }
    catch( profugus::assertion const & a )
    {
        error_message = std::string( a.what() );
        caught = true;
    }
    catch( ... )
    {
        caught = false;
    }

    EXPECT_TRUE(caught);

    // Make sure we can extract the error message.
    std::string const compare_value(
        "Assertion: hello1, failed in myfile:42" );
    EXPECT_EQ(0, error_message.compare( compare_value ));
}

//---------------------------------------------------------------------------//

TEST(Catch, string_literal)
{
    bool caught = false;
    try
    {
        throw "hello";
    }
    catch( profugus::assertion const & a )
    {
        caught = false;
    }
    catch( const char* msg )
    {
        caught = true;
    }
    catch( ... )
    {
        caught = false;
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//

TEST(DBC, toss_cookies)
{
    bool caught = false;
    try
    {
        std::string const msg("testing toss_cookies()");
        std::string const file("DummyFile.ext");
        int const line( 55 );
        profugus::toss_cookies( msg.c_str(), file.c_str(), line );
        throw "Bogus!";
    }
    catch( profugus::assertion const & a )
    {
        caught = true;
    }
    catch( ... )
    {
        caught = false;
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//

TEST(DBC, Require)
{
    bool caught = false;

    try
    {
        REQUIRE( 0 );
        throw "Bogus!";
    }
    catch( profugus::assertion const & a )
    {
#if UTILS_DBC & 1
        caught = true;
#else
        caught = false;
#endif
    }
    catch(...)
    {
#if UTILS_DBC & 1
        caught = false;
#else
        caught = true;
#endif
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//

TEST(DBC, Check)
{
    bool caught = false;

    try
    {
        CHECK( 0 );
        throw "Bogus!";
    }
    catch( profugus::assertion const & a )
    {
#if UTILS_DBC & 2
        caught = true;
#else
        caught = false;
#endif
    }
    catch(...)
    {
#if UTILS_DBC & 2
        caught = false;
#else
        caught = true;
#endif
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//

TEST(DBC, Ensure)
{
    bool caught = false;

    try
    {
        ENSURE( 0 );
        throw "Bogus!";
    }
    catch( profugus::assertion const & a )
    {
#if UTILS_DBC & 4
        caught = true;
#else
        caught = false;
#endif
    }
    catch(...)
    {
#if UTILS_DBC & 4
        caught = false;
#else
        caught = true;
#endif
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//

TEST(DBC, Insist)
{
    bool caught = false;

    try
    {
        Insist( 0 , "This is a bad thing.");
        throw "Bogus!";
    }
    catch( profugus::assertion const & a )
    {
        caught = true;
    }
    catch(...)
    {
        caught = false;
    }

    EXPECT_TRUE(caught);
}

//---------------------------------------------------------------------------//
//                 end of tstDBC.cc
//---------------------------------------------------------------------------//
