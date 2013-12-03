//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Scalar_Unit_Test.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:54:57 2008
 * \brief  Scalar_Unit_Test member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include "DBC.hh"
#include "Scalar_Unit_Test.hh"

namespace nemesis
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for Scalar_Unit_Test
 *
 * \arg argc The number of command line arguments
 * \arg argv A list of strings containg the command line arguments
 * \arg release_ A function pointer to this package's release function.
 * \arg out_ A user specified iostream that defaults to std::cout.
 * \exception nemesis::assertion An exception with the message "Success" will
 * be thrown if \c --version is found in the argument list.
 *
 * The constructor initializes the base class UnitTest by setting numPasses
 * and numFails to zero.  It also prints a message that declares this to be a
 * scalar unit test and provides the unit test name.
 */
Scalar_Unit_Test::Scalar_Unit_Test(int &argc,
                                   char **&argv,
                                    string_fp_void release_,
                                    std::ostream & out_ )
    : Unit_Test( argc, argv, release_, out_ )
{
    using std::endl;
    using std::string;

    Require( argc > 0 );
    Require( release != NULL );

    // header

    out << "\n============================================="
        << "\n=== Scalar Unit Test: " << testName
        << "\n=============================================\n" << endl;

    // version tag

    out << testName << ": version " << release() << "\n" << endl;

    // exit if command line contains "--version"

    for( int arg = 1; arg < argc; arg++ )
        if( string( argv[arg] ) == "--version" )
            throw nemesis::assertion( string( "Success" ) );

    Ensure( numPasses == 0 );
    Ensure( numFails  == 0 );
    Ensure( testName.length() > 0 );

    return;
}

} // end namespace nemesis

//---------------------------------------------------------------------------//
//                 end of Scalar_Unit_Test.cc
//---------------------------------------------------------------------------//
