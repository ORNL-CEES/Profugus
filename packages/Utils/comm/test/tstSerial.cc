//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstSerial.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:37:03 2008
 * \brief  Serial comm test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "release/Release.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstScalar(Parallel_Unit_Test &ut)
{

// Skip the tests if code not configured with the option --with-comm=scalar.

#ifndef COMM_SCALAR

    if( nemesis::isScalar() )
        ut.failure("Incorrectly identified process as scalar.");
    else
        ut.passes("Correctly identified process as parallel.");

#else

    // Check the isScalar function.
    if( nemesis::isScalar() )
        ut.passes("Correctly identified process as scalar.");
    else
        ut.failure("Incorrectly identified process as parallel.");

    // For --with-comm=scalar, probe(int,int,int) always returns false.

    int int3(99);
    bool const probeResult( nemesis::probe(0,0,int3) );
    if( probeResult )
    {
        ut.failure( "For --with-comm=scalar, probe(int,int,int) returned true." );
    }
    else
    {
        ut.passes( "For --with-comm=scalar, probe(int,int,int) returned false." );
    }

    // Test broadcast function

    int result(0);
    result = nemesis::broadcast( &int3, 1, 0 );
    if( result == nemesis::COMM_SUCCESS )
    {
        ut.passes("For --with-comm=scalar, broadcast() returned COMM_SUCCCESS.");
    }
    else
    {
        ut.passes("For --with-comm=scalar, broadcast() did not return COMM_SUCCCESS.");
    }

    // BLOCKING send/receive not allowed
    bool caught = false;
    try
    {
        result = nemesis::send( &int3, 1, 0, 0 );
    }
    catch (const nemesis::assertion &a)
    {
        caught = true;
    }
    if (caught)
    {
        ut.passes("For --with-comm=scalar, send() disallowed.");
    }
    else
    {
        ITFAILS;
    }

    caught = false;
    try
    {
        result = nemesis::receive( &int3, 1, 0, 0 );
    }
    catch (const nemesis::assertion &a)
    {
        caught = true;
    }
    if (caught)
    {
        ut.passes("For --with-comm=scalar, receive() disallowed.");
    }
    else
    {
        ITFAILS;
    }

#endif
}

//---------------------------------------------------------------------------//

void test_self_comm(Parallel_Unit_Test &ut)
{
#ifdef COMM_SCALAR

    // make a receive buffer
    int size = 0;

    // post a receive
    nemesis::Request r = nemesis::receive_async(&size, 1, 0, 101);

    // send the size
    int send_size = 6;
    nemesis::Request s = nemesis::send_async(&send_size, 1, 0, 101);

    UNIT_TEST(r.inuse());
    UNIT_TEST(s.inuse());

    // wait on the message
    s.wait();
    r.wait();

    UNIT_TEST(!r.inuse());
    UNIT_TEST(!s.inuse());

    UNIT_TEST(size == 6);

    // post receives
    vector<double>   receive_buffer(size, 0.0);
    nemesis::Request r2;
    nemesis::receive_async(r2, &receive_buffer[0], size, 0, 201);

    for (int i = 0; i < size; ++i)
    {
        UNIT_TEST(receive_buffer[i] == 0.0);
    }

    // send sized buffer
    vector<double> data(send_size, 0.1);
    data[3] = 0.2;
    data[5] = 0.9;
    nemesis::Request s2;
    nemesis::send_async(s2, &data[0], send_size, 0, 201);

    UNIT_TEST(r2.inuse());
    UNIT_TEST(s2.inuse());

    // wait on the message
    s2.wait();
    r2.wait();

    UNIT_TEST(!r2.inuse());
    UNIT_TEST(!s2.inuse());

    UNIT_TEST(receive_buffer[0] == 0.1);
    UNIT_TEST(receive_buffer[1] == 0.1);
    UNIT_TEST(receive_buffer[2] == 0.1);
    UNIT_TEST(receive_buffer[3] == 0.2);
    UNIT_TEST(receive_buffer[4] == 0.1);
    UNIT_TEST(receive_buffer[5] == 0.9);

    if (ut.numFails == 0)
        ut.passes("Non-blocking comm-to-self works correctly.");
#endif
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    nemesis::Parallel_Unit_Test ut(argc, argv, nemesis::release::short_version);

    node  = nemesis::node();
    nodes = nemesis::nodes();

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        tstScalar(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        test_self_comm(ut);
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
        std::cout << "ERROR: While testing tstSerial, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSerial, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstSerial.cc
//---------------------------------------------------------------------------//
