//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstDBC.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 12:44:57 2008
 * \brief  DBC tests.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "Release.hh"
#include "../DBC.hh"

int num_failed = 0;
int num_passed = 0;

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
//The way this test article works is that each of the DBC macros are tested in
//a seperate function.  A false condition is asserted using each macro, and
//after this follows a throw.  Two catch clauses are available, one to catch
//an assertion object, and one to catch anything else.  By comparing the
//exception that is actually caught with the one that should be caught given
//the DBC setting in force, we can determine whether each test passes or
//fails.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Make sure we can differentiate betweeen a std::runtime_error and a
// nemesis::assertion.
//---------------------------------------------------------------------------//

static void t1()
{
    std::cout << "t1 Test: ";
    try
    {
        throw std::runtime_error( "hello1" );
    }
    catch( nemesis::assertion const & a )
    {
        std::cout << a.what() << std::endl;
        std::cout << "failed" << std::endl;
        ++num_failed;
    }
    catch( ... )
    {
        std::cout << "passed" << std::endl;
        ++num_passed;
    }
    return;
}

//---------------------------------------------------------------------------//
// Make sure we can catch a nemesis::assertion and extract the error
// message.
//---------------------------------------------------------------------------//

static void t2()
{
    std::cout << "t2-a Test: ";
    std::string error_message;
    try
    {
        throw nemesis::assertion( "hello1", "myfile", 42 );
    }
    catch( nemesis::assertion const & a )
    {
        std::cout << "passed" << std::endl;
        error_message = std::string( a.what() );
        ++num_passed;
    }
    catch( ... )
    {
        std::cout << "failed" << std::endl;
        ++num_failed;
    }

    // Make sure we can extract the error message.

    std::cout << "t2-b Test: ";
    std::string const compare_value(
        "Assertion: hello1, failed in myfile:42" );
    if ( error_message.compare( compare_value ) == 0 )
    {
        std::cout << "passed" << std::endl;
        ++num_passed;
    }
    else
    {
        std::cout << "failed" << std::endl;
        ++num_failed;
    }
    return;
}

//---------------------------------------------------------------------------//
// Test throwing and catching of a literal
//
//lint -e1775  do not warn about "const char*" not being a declared exception
//             type.
//lint -e1752  do not warn about catching "const char*" instead of
//             catching a reference to an exception (see "More Effective C++"
//             for details about catching exceptions)
//---------------------------------------------------------------------------//

static void t3()
{
    std::cout << "t3 Test: ";
    try
    {
        throw "hello";
    }
    catch( nemesis::assertion const & a )
    {
        std::cout << a.what() << std::endl;
        std::cout << "failed" << std::endl;
        ++num_failed;
    }
    catch( const char* msg )
    {
        std::cout << "passed   "
                  << "msg = " << msg << std::endl;
        ++num_passed;
    }
    catch( ... )
    {
        std::cout << "failed" << std::endl;
        ++num_failed;
    }
    return;
}

//---------------------------------------------------------------------------//
// Check the toss_cookies function.
// This function builds an error message and throws an exception.
//---------------------------------------------------------------------------//

static void ttoss_cookies()
{
    {
        std::cout << "ttoss_cookies Test: ";
        try
        {
            std::string const msg("testing toss_cookies()");
            std::string const file("DummyFile.ext");
            int const line( 55 );
            nemesis::toss_cookies( msg.c_str(), file.c_str(), line );
            throw "Bogus!";
        }
        catch( nemesis::assertion const & a )
        {
            std::cout << "passed" << std::endl;
            ++num_passed;
        }
        catch( ... )
        {
            std::cout << "failed" << std::endl;
            ++num_failed;
        }
    }
    return;
}

//---------------------------------------------------------------------------//
// Check the operation of the Require() macro.
//---------------------------------------------------------------------------//

static void trequire()
{
    std::cout << "t-Require Test: ";
    try {
        Require( 0 );
        throw "Bogus!";
    }
    catch( nemesis::assertion const & a )
    {
#if NEMESIS_DBC & 1
        std::cout << "passed" << std::endl;
        ++num_passed;

        std::cout << "t-Require message value Test: ";
        {
            std::string msg( a.what() );
            std::string expected_value( "Assertion: 0, failed in" );
            string::size_type idx = msg.find( expected_value );
            if( idx != string::npos )
            {
                cout << "passed" << std::endl;
                ++num_passed;
            }
            else
            {
                cout << "failed" << std::endl;
                ++num_failed;
            }
        }

#else
        std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
        ++num_failed;
#endif
    }
    catch(...)
    {
#if NEMESIS_DBC & 1
        std::cout << "failed" << std::endl;
        ++num_failed;
#else
        std::cout << "passed" << std::endl;
        ++num_passed;
#endif
    }
    return;
}

//---------------------------------------------------------------------------//
// Check the operation of the Check() macro.
//---------------------------------------------------------------------------//

static void tcheck()
{
    std::cout << "t-Check Test: ";
    try {
        Check( false );
        throw std::runtime_error( std::string( "tstAssert: t2()" ) );
    }
    catch( nemesis::assertion const & a )
    {
#if NEMESIS_DBC & 2
        std::cout << "passed" << std::endl;
        ++num_passed;

        std::cout << "t-Check message value Test: ";
        {
            std::string msg( a.what() );
            std::string expected_value( "Assertion: false, failed in" );
            string::size_type idx = msg.find( expected_value );
            if( idx != string::npos )
            {
                cout << "passed" << std::endl;
                ++num_passed;
            }
            else
            {
                cout << "failed" << std::endl;
                ++num_failed;
            }
        }
#else
        std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
        std::string msg( a.what() );
        ++num_failed;
#endif
    }
    catch(...)
    {
#if NEMESIS_DBC & 2
        std::cout << "failed\n";
        ++num_failed;
#else
        std::cout << "passed\n";
        ++num_passed;
#endif
    }
    return;
}

//---------------------------------------------------------------------------//
// Check the operation of the Ensure() macro.
//---------------------------------------------------------------------------//

static void tensure()
{
    std::cout << "t-Ensure Test: ";
    try {
        Ensure(0);
        throw "Bogus!";
    }
    catch( nemesis::assertion const & a )
    {
#if NEMESIS_DBC & 4
        std::cout << "passed" << std::endl;
        ++num_passed;

        std::cout << "t-Ensure message value Test: ";
        {
            std::string msg( a.what() );
            std::string expected_value( "Assertion: 0, failed in" );
            string::size_type idx = msg.find( expected_value );
            if( idx != string::npos )
            {
                cout << "passed" << std::endl;
                ++num_passed;
            }
            else
            {
                cout << "failed" << std::endl;
                ++num_failed;
            }
        }

#else
        std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
#endif
    }
    catch(...)
    {
#if NEMESIS_DBC & 4
        std::cout << "failed\n";
        ++num_failed;
#else
        std::cout << "passed\n";
        ++num_passed;
#endif
    }
    return;
}

//---------------------------------------------------------------------------//

static void tremember()
{
    //lint -e774  do not warn about if tests always evaluating to False.  The
    //            #if confuses flexelint here.

    std::cout << "t-Remember Test: ";

    int x = 0;
    Remember(x = 5);
#if NEMESIS_DBC & 4
    if (x != 5)
    {
        std::cout << "failed" << std::endl;
        ++num_failed;
    }
    else
    {
        std::cout << "passed" << std::endl;
        ++num_passed;
    }
#else
    if (x != 0)
    {
        std::cout << "failed" << std::endl;
        ++num_failed;
    }
    else
    {
        std::cout << "passed" << std::endl;
        ++num_passed;
    }
#endif
    return;}

//---------------------------------------------------------------------------//
// Check the operation of the Assert() macro, which works like Check().
//---------------------------------------------------------------------------//

static void tassert()
{
    std::cout << "t-Assert Test: ";
    try {
        Assert(0);
        throw "Bogus!";
    }
    catch( nemesis::assertion const & a )
    {
#if NEMESIS_DBC & 2
        std::cout << "passed" << std::endl;
        ++num_passed;

        std::cout << "t-Assert message value Test: ";
        {
            std::string msg( a.what() );
            std::string expected_value( "Assertion: 0, failed in" );
            string::size_type idx = msg.find( expected_value );
            if( idx != string::npos )
            {
                cout << "passed" << std::endl;
                ++num_passed;
            }
            else
            {
                cout << "failed" << std::endl;
                ++num_failed;
            }
        }
#else
        std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
        ++num_failed;
#endif
    }
    catch(...)
    {
#if NEMESIS_DBC & 2
        std::cout << "failed\n";
        ++num_failed;
#else
        std::cout << "passed\n";
        ++num_passed;
#endif
    }
    return;
}

//---------------------------------------------------------------------------//
// Basic test of the Insist() macro.
//---------------------------------------------------------------------------//

static void tinsist()
{
    //lint -e506  Do not warn about constant value boolean in the Insist
    //            test.
    {
        std::cout << "t-Insist Test: ";
        std::string insist_message( "You must be kidding!" );
        try {
            Insist( 0, insist_message );
            throw "Bogus!";
        }
        catch( nemesis::assertion const & a )
        {
            std::cout << "passed" << std::endl;
            ++num_passed;

            std::cout << "t-Insist message value Test: ";
            {
                bool passed( true );
                std::string msg( a.what() );
                std::string expected_value( "You must be kidding!" );
                string::size_type idx( msg.find( expected_value ) );
                if( idx == string::npos ) passed=false;
                idx = msg.find( insist_message );
                if( idx == string::npos ) passed=false;
                if( passed )
                {
                    cout << "passed" << std::endl;
                    ++num_passed;
                }
                else
                {
                    cout << "failed" << std::endl;
                    ++num_failed;
                }
            }
        }
        catch(...)
        {
            std::cout << "failed: wrong exception caught" << std::endl;
            ++num_failed;
        }
    }

    return;
}

//---------------------------------------------------------------------------//
// Basic test of the NotImplemented() macro.
//---------------------------------------------------------------------------//

static void tnotimpl()
{
    std::cout << "t-NotImplemented Test: ";
    std::string nimplmsg( "frozzing" );
    try {
        Not_Implemented(nimplmsg.c_str());
        throw "Bogus!";
    }
    catch( nemesis::not_implemented_error const & a )
    {
        std::cout << "passed" << std::endl;
        ++num_passed;

        std::cout << "t-NotImplemented message value Test: ";
        {
            bool passed( true );
            std::string msg( a.what() );
            string::size_type idx( msg.find( nimplmsg ) );
            if( idx == string::npos ) passed=false;
            if( passed )
            {
                cout << "passed" << std::endl;
                ++num_passed;
            }
            else
            {
                ++num_failed;
            }
        }
    }
    catch(...)
    {
        std::cout << "failed: wrong exception caught" << std::endl;
        ++num_failed;
    }

    return;
}

//---------------------------------------------------------------------------//
// Basic test of the Validate() macro.
//---------------------------------------------------------------------------//

static void tvalidate()
{
    std::cout << "t-Validate Test: ";
    try {
        Validate(2 + 2 > 4,
                "Math works like " << 2 << " and " << 4);
        throw "Bogus!";
    }
    catch( nemesis::validation_error const & a )
    {
        std::cout << "passed" << std::endl;
        ++num_passed;

        std::cout << "t-Validate message value Test: ";
        {
            bool passed( true );
            std::string msg( a.what() );
            cout << "\n\n\n" << msg << endl;
            std::string expected("Math works like 2 and 4");
            string::size_type idx( msg.find( expected ) );
            if( idx == string::npos ) passed=false;
            if( passed )
            {
                cout << "passed" << std::endl;
                ++num_passed;
            }
            else
            {
                cout << "failed" << std::endl;
                ++num_failed;
            }
        }
    }
    catch(...)
    {
        std::cout << "failed: wrong exception caught" << std::endl;
        ++num_failed;
    }

    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if( string( argv[arg] ).find( "--version" ) == 0 )
        {
            cout << argv[0] << ": version " << nemesis::release()
                 << endl;
            return 0;
        }

    // >>> UNIT TESTS

    // Test basic throw and catch functionality.
    t1();
    t2();
    t3();

    // Test mechanics of Assert funtions.
    ttoss_cookies();

    // Test Design-by-Constract macros.
    trequire();
    tcheck();
    tensure();
    tremember();
    tassert();
    tinsist();
    tnotimpl();
    tvalidate();

    cout << endl;
    cout << "\n*********************************************\n";
    if (num_passed > 0 && num_failed == 0)
    {
        cout << "**** tstDBC Test: PASSED" << endl;
    }
    else
    {
        cout << "**** tstDBC Test: FAILED" << endl;
    }
    cout << "*********************************************\n";
}

//---------------------------------------------------------------------------//
//                        end of tstDBC.cc
//---------------------------------------------------------------------------//
