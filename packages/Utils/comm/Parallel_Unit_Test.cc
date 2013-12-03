//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Parallel_Unit_Test.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 15:44:20 2008
 * \brief  Parallel_Unit_Test class definition.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include "SpinLock.hh"
#include "global.hh"
#include "Parallel_Unit_Test.hh"

namespace nemesis
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for Parallel_Unit_Test
 *
 * \arg argc The number of command line arguments
 * \arg argv A list of strings containg the command line arguments
 * \arg release_ A function pointer to this package's release function.
 * \arg out_ A user specified iostream that defaults to std::cout.
 *
 * \exception nemesis::assertion An exception with the message "Success" will
 * be thrown if \c --version is found in the argument list.
 *
 * The constructor initializes the prallel communicator (MPI) and then
 * initializes the base class UnitTest by setting numPasses and numFails to
 * zero.  It also prints a message that declares this to be a scalar unit test
 * and provides the unit test name.
 */
Parallel_Unit_Test::Parallel_Unit_Test(int              &argc,
                                       char           **&argv,
                                       string_fp_void    release_,
                                       std::ostream     &out_ )
    : Unit_Test( argc, argv, release_, out_ )
{
    using std::string;

    initialize( argc, argv );

    Require( argc > 0 );
    Require( release != NULL );

    // header

    if( node() == 0 )
        out << "\n============================================="
            << "\n=== Parallel Unit Test: " << testName
            << "\n=== Number of Processors: " << nodes()
            << "\n=============================================\n"
            << std::endl;

    // version tag

    if( node() == 0 )
        out << testName << ": version " << release() << "\n" << std::endl;

    // exit if command line contains "--version"

    for( int arg = 1; arg < argc; arg++ )
        if( string( argv[arg] ) == "--version" )
        {
            finalize();
            throw assertion( string( "Success" ) );
        }
    Ensure( numPasses == 0 );
    Ensure( numFails  == 0 );
    Ensure( testName.length() > 0 );

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * The destructor provides a final status report before it calls MPI_Finalize
 * and exits.
 */
Parallel_Unit_Test::~Parallel_Unit_Test()
{
    global_barrier();
    if( node() == 0 )
        out << resultMessage() << std::endl;
    global_barrier();
    finalize();
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print a summary of the pass/fail status of Parallel_Unit_Test.
 */
void Parallel_Unit_Test::status()
{
    {
        // Provide some space before the report -- but keep all the processors
        // in sync.
        HTSyncSpinLock slock;
        if( node() == 0 )
            out << std::endl;
    }
    {
        HTSyncSpinLock slock;
        out << "Done testing " << testName << " on node "
            << node() << "." << std::endl;
    }
    return;
}

} // end namespace nemesis

//---------------------------------------------------------------------------//
//                 end of Parallel_Unit_Test.cc
//---------------------------------------------------------------------------//
