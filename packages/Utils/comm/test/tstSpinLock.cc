//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstSpinLock.cc
 * \author Thomas M. Evans
 * \date   Tue Jun 07 13:17:34 2011
 * \brief  SpinLock unit test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "release/Release.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"
#include "../SpinLock.hh"

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

void ht_test(Parallel_Unit_Test &ut)
{
    // make file
    if (node == 0)
    {
        ofstream out("ht.out");
    }

    // make data on each node
    vector<int> x(5, node);
    for (int i = 1; i < 5; ++i)
        x[i] = x[i - 1] + 1;

    {
        nemesis::HTSyncSpinLock l;

        cout << ">>> Spin-lock writing output on node " << node << endl;

        ofstream out("ht.out", ios::app);
        for (int i = 0; i < 5; ++i)
            out << x[i] << " ";
        out << endl;
    }

    // check output
    if (node == 0)
    {
        ifstream in("ht.out");
        for (int n = 0; n < nodes; ++n)
        {
            for (int i = 0; i < 5; ++i)
            {
                int ref;
                in >> ref;
                UNIT_TEST(ref == i + n);
            }
        }

        if (ut.numFails == 0)
        {
            ut.passes("Head-tail spin-locking appears to operate correctly");
        }
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
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        ht_test(ut);
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
        std::cout << "ERROR: While testing tstSpinLock, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSpinLock, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstSpinLock.cc
//---------------------------------------------------------------------------//
