//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstAbort.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:34:00 2008
 * \brief  Abort test.
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

void abort_test(Parallel_Unit_Test &ut)
{
    cout << "Entering abort on proc " << nemesis::node() << endl;

    nemesis::global_barrier();

    // only abort from processor 0 for nice output
    if (nemesis::node() == 0)
    {
        cout << "Aborting from processor 0" << endl;
        nemesis::abort();
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
        int gpass = 1;
        int gfail = 0;

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstAbort, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstAbort, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstAbort.cc
//---------------------------------------------------------------------------//
