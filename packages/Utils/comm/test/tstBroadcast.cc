//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstBroadcast.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:34:18 2008
 * \brief  Broadcast test.
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

using namespace std;
using nemesis::soft_equiv;
using nemesis::Parallel_Unit_Test;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_simple(Parallel_Unit_Test &ut)
{
    using nemesis::broadcast;

    char   c = 0;
    int    i = 0;
    long   l = 0;
    float  f = 0;
    double d = 0;
    vector<double> vref(10,3.1415);
    vector<double> v(10,0.0);

    // assign on node 0
    if (node == 0)
    {
        c = 'A';
        i = 1;
        l = 1000;
        f = 1.5;
        d = 2.5;
        v = vref;
    }

    // send out data, using node 0 as root
    broadcast(&c, 1, 0);
    broadcast(&i, 1, 0);
    broadcast(&l, 1, 0);
    broadcast(&f, 1, 0);
    broadcast(&d, 1, 0);
    broadcast(v.begin(),v.end(),v.begin());

    // check values
    if (c != 'A')             ITFAILS;
    if (i != 1)               ITFAILS;
    if (l != 1000)            ITFAILS;
    if (!soft_equiv(f, 1.5f)) ITFAILS;
    if (!soft_equiv(d, 2.5))  ITFAILS;
    if (!soft_equiv(v.begin(),v.end(),vref.begin(),vref.end()))  ITFAILS;

    nemesis::global_barrier();

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "test_simple() ok on " << nemesis::node();
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//
//  By adjusting the parameters below, this test will overflow the MPI memory
//  buffers.  Read the comments below if you'd like to do this.

void test_loop(Parallel_Unit_Test &ut)
{
    using nemesis::broadcast;

    // >>> kmax controls how much data is broadcast.  If kmax is too big
    // >>> (like 10000000), shmem will fail.
    const int kmax = 10;

    if (nemesis::node() == 0) // host proc
    {
        // send out the values on host
        for ( int k = 0; k < kmax; ++k )
        {
            Insist(! broadcast(&k, 1, 0), "MPI Error");
            double foo = k + 0.5;
            Insist(! broadcast(&foo, 1, 0), "MPI Error");
        }
    }
    else // all other procs
    {
        // >>> Use sleep() if you want the host processor to fill up the
        // >>> buffers.  We comment out the sleep() command here because it's
        // >>> not supported on all all platforms.

        // sleep(10);

        int kk;
        double foofoo;
        for ( int k = 0; k < kmax; ++k )
        {
            kk = -1;
            foofoo = -2.0;
            Insist(! broadcast(&kk, 1, 0), "MPI Error");
            if ( kk != k ) ITFAILS;
            Insist(! broadcast(&foofoo, 1, 0), "MPI Error");
            if ( foofoo != k + 0.5 ) ITFAILS;
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "test_loop() ok on " << nemesis::node();
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, nemesis::release::short_version);

    node  = nemesis::node();
    nodes = nemesis::nodes();

    try
    {
        int gpass = 0;
        int gfail = 0;

        // >>> UNIT TESTS
        test_simple(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        test_loop(ut);
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
        std::cout << "ERROR: While testing tstBroadcast, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstBroadcast, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstBroadcast.cc
//---------------------------------------------------------------------------//
