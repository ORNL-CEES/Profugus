//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstGather.cc
 * \author Thomas M. Evans
 * \date   Mon Jul 15 12:23:15 2013
 * \brief  Gather comm tests.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "release/Release.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"

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

void all_gather_test(Parallel_Unit_Test &ut)
{
    const int num_els = 2;
    vector<int> inp(num_els);
    vector<int> outp(num_els * nodes, -1);

    // Set data to send
    inp[0] = 10 * node;
    inp[1] = 10 * node + 1;

    // nemesis::global_barrier();
    // cout << "Input on node " << node << ": ";
    // copy(inp.begin(), inp.end(), ostream_iterator<int>(cout, ", "));
    // cout << endl;

    // Do communication
    nemesis::all_gather(&inp[0], &outp[0], num_els);

    // cout << "Output on node " << node << ": ";
    // copy(outp.begin(), outp.end(), ostream_iterator<int>(cout, ", "));
    // cout << endl;
    // nemesis::global_barrier();

    // Check received data
    for (int n = 0; n < nodes; ++n)
    {
        UNIT_TEST(outp[2*n + 0] == 10 * n + 0);
        UNIT_TEST(outp[2*n + 1] == 10 * n + 1);
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "All_gather ok on " << node << " of " << nodes;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

void gather_test(Parallel_Unit_Test &ut)
{
    const int num_els = 2;
    vector<int>  inp(num_els);
    vector<int>  outp;
    int         *outp_ptr = 0;

    if (node == 0)
    {
        outp.resize(num_els * nodes, -1);
        outp_ptr = &outp[0];
    }

    // Set data to send
    inp[0] = 10 * node;
    inp[1] = 10 * node + 1;

    // gather to node 0
    nemesis::gather(&inp[0], outp_ptr, num_els, 0);

    if (node == 0)
    {
        for (int n = 0; n < nodes; ++n)
        {
            UNIT_TEST(outp[2*n + 0] == 10 * n + 0);
            UNIT_TEST(outp[2*n + 1] == 10 * n + 1);
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Gather ok on " << node << " of " << nodes;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, denovo::release);

    node  = nemesis::node();
    nodes = nemesis::nodes();

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        all_gather_test(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        gather_test(ut);
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
        std::cout << "ERROR: While testing tstGather, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstGather, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstGather.cc
//---------------------------------------------------------------------------//
