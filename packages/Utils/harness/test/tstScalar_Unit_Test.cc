//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstScalar_Unit_Test.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 12:58:24 2008
 * \brief  Unit test for testing harness.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <map>
#include <cstdlib>

#include "Release.hh"
#include "../Scalar_Unit_Test.hh"

using namespace std;
using namespace nemesis;

// Provide old style call to pass/fail macros.  Use object name unitTest for
// this unit test.
#define PASSMSG(a) unitTest.passes(a)
#define ITFAILS    unitTest.failure(__LINE__);
#define FAILURE    unitTest.failure(__LINE__, __FILE__);
#define FAILMSG(a) unitTest.failure(a);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstOne( Unit_Test &unitTest )
{
    unitTest.passes("Looks like the passes member function is working.");
    PASSMSG("Looks like the PASSMSG macro is working as a member function.");

    return;
}

//---------------------------------------------------------------------------//

void tstTwo( Unit_Test &unitTest )
{
    unitTest.failure("Looks like the failure member function is working.");
    FAILMSG("Looks like the FAILMSG macro is working.");
    ITFAILS;
    FAILURE;

    // Kill report of failures
    unitTest.reset();

    // We need at least one pass.
    PASSMSG("Done with tstTwo.");
    return;
}

//---------------------------------------------------------------------------//

void tstTwoCheck( Unit_Test &unitTest, std::ostringstream & msg )
{
    bool verbose(true);
    std::map<string,unsigned> word_list(
        Unit_Test::get_word_count( msg, verbose ) );

    // Check the list of occurances against the expected values
    if( word_list[ string("Test") ] == 6 )
        unitTest.passes("Found 6 occurances of \"Test\"");
    else
        unitTest.failure("Did not find expected number of occurances of \"Test\"");

    if( word_list[ string("failed") ] != 4 )
        unitTest.failure("Did not find 4 occurances of failure.");
    if( word_list[ string("FAILMSG") ] != 1 )
        unitTest.failure("Found 1 occurance of \"FAILMSG\"");
    if( word_list[ string("failure") ] != 1 )
        unitTest.failure("Found 1 occurance of \"failure\"");

    if( word_list[ string("macro") ] == 1 )
        unitTest.passes("Found 1 occurance of \"macro\"");
    else
        unitTest.failure("Did not find expected number of occurances of \"macro\"");

    if( word_list[ string("working") ] == 2 )
        unitTest.passes("Found 2 occurances of \"working\"");
    else
        unitTest.failure("Did not find expected number of occurances of \"working\"");

    return;
}

//---------------------------------------------------------------------------//

void tstVersion( Unit_Test & ut, int & argc, char **& argv )
{
    // build the command that contains "--version"
    string cmd;
    for( int ic=0; ic<argc; ++ic )
        cmd += " " + string( argv[0] );
    cmd += " --version";

    system( cmd.c_str() );
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        // >>> UNIT TESTS

        // Test ctor for Scalar_Unit_Test (also tests Unit_Test ctor and member
        // function setTestName).
        Scalar_Unit_Test ut( argc, argv, release );
        tstOne(ut);

        // Silent version.
        std::ostringstream messages;
        Scalar_Unit_Test sut( argc, argv, release, messages );
        tstTwo(sut);

        tstTwoCheck( ut, messages );
        if( argc == 1 )
        {
            // Test --version option.
            tstVersion( ut, argc, argv );
        }
    }
    catch(assertion &err)
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (std::exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << "An unknown exception was thrown" << endl;
        return 1;
    }

    return 0;
}

//---------------------------------------------------------------------------//
//                        end of tstScalar_Unit_Test.cc
//---------------------------------------------------------------------------//
