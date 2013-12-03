//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstTime.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:37:16 2008
 * \brief  Timer test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "release/Release.hh"
#include "harness/Soft_Equivalence.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"
#include "../Timer.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void wall_clock_test(Parallel_Unit_Test &ut)
{
    using std::endl;
    using std::cout;
    using std::ostringstream;

    using nemesis::wall_clock_time;
    using nemesis::wall_clock_resolution;
    using nemesis::Timer;

    double const wcr( nemesis::wall_clock_resolution() );
    if( wcr > 0.0 && wcr <= 100.0)
    {
        ostringstream msg;
        msg << "The timer has a wall clock resoution of "
            << wcr << " ticks." << endl;
        ut.passes(msg.str());
    }
    else
    {
        ostringstream msg;
        msg << "The timer does not appear to have a reasonable resolution."
            << " nemesis::wall_clock_resolution() = " << wcr << " ticks."
            << endl;
        ut.failure(msg.str());
    }

    Timer t;
    double begin = nemesis::wall_clock_time();
    t.start();

    double const prec( 1.5*t.posix_err() );


    for( int i = 0; i < 200000000; i++ )
    { /* empty */
    }

    double end = nemesis::wall_clock_time();
    t.stop();

    double const error( t.wall_clock() - (end-begin) );
    if( std::fabs(error) <= prec )
    {
        ut.passes("wall_clock() value looks ok.");
    }
    else
    {
        ostringstream msg;
        msg << "t.wall_clock() value does not match the expected value."
            << "\n\tend            = " << end
            << "\n\tbegin          = " << begin
            << "\n\tend-begin      = " << end - begin
            << "\n\tt.wall_clock() = " << t.wall_clock()
            << "\n\tprec           = " << prec << endl;
        ut.failure(msg.str());
    }

    //---------------------------------------------------------------------//
    // Ensure that system + user <= wall
    //
    // Due to round off errors, the wall clock time might be less than the
    // system + user time.  But this difference should never exceed
    // t.posix_err().
    //---------------------------------------------------------------------//

    double const deltaWallTime( t.wall_clock() - (
                                    t.system_cpu() + t.user_cpu() ) );

    if( deltaWallTime > 0.0 || std::fabs(deltaWallTime) <= prec )
    {
        ostringstream msg;
        msg << "The sum of cpu and user time is less than or equal to the\n\t"
            << "reported wall clock time (within error bars = " << prec
            << " secs.)." << endl;
        ut.passes(msg.str());
    }
    else
    {
        ostringstream msg;
        msg << "The sum of cpu and user time exceeds the reported wall "
            << "clock time.  Here are the details:"
            << "\n\tposix_error() = " << prec << " sec."
            << "\n\tdeltaWallTime = " << deltaWallTime  << " sec."
             << "\n\tSystem time   = " << t.system_cpu() << " sec."
             << "\n\tUser time     = " << t.user_cpu()   << " sec."
             << "\n\tWall time     = " << t.wall_clock() << " sec."
            << endl;
        ut.failure(msg.str());
    }

    //------------------------------------------------------//
    // Demonstrate print functions:
    //------------------------------------------------------//

    cout << "\nDemonstration of the print() member function via the\n"
         << "\toperator<<(ostream&,Timer&) overloaded operator.\n"
         << endl;

    cout << "Timer = " << t << endl;

    //------------------------------------------------------//
    // Do a second timing:
    //------------------------------------------------------//

    cout << "\nCreate a Timer Report after two timing cycles:\n"
         << endl;

    t.start();
    for( int i = 0; i < 200000000; i++ )
    { /* empty */
    }
    t.stop();

    t.print( cout, 6 );

    //------------------------------------------------------//
    // Check the number of intervals
    //------------------------------------------------------//

    int const expectedNumberOfIntervals(2);
    if( t.intervals() == expectedNumberOfIntervals )
        ut.passes("Found the expected number of intervals.");
    else
        ut.failure("Did not find the expected number of intervals.");

    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, nemesis::release::short_version);

    node  = nemesis::node();
    nodes = nemesis::nodes();

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        wall_clock_test(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstTime, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTime, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstTime.cc
//---------------------------------------------------------------------------//
