//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstParallel_Utils.cc
 * \author Thomas M. Evans
 * \date   Fri Jul 13 15:26:57 2007
 * \brief  Parallel_Utils test.
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
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
#include "../Parallel_Utils.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);

int global_passes = 0;
int global_fails  = 0;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void equivalence_check(Parallel_Unit_Test &ut)
{
    // define some stuff on each pe
    int x    = 12;
    double y = 0.001;

    // everybody should be the same
    if (!nemesis::check_global_equiv(x)) ITFAILS;
    if (!nemesis::check_global_equiv(y)) ITFAILS;

    if (nodes > 1)
    {
        if (node == 1)
        {
            x = 13;
            y = 0.00101;
        }

        // the test will always fail on processor 0, and it will fail on
        // processor 1 when nodes > 2

        bool result_x = nemesis::check_global_equiv(x);
        bool result_y = nemesis::check_global_equiv(y);

        if (node == 0)
        {
            if (result_x) ITFAILS;
            if (result_y) ITFAILS;
        }

        if (nodes > 2)
        {
            if (node == 1)
            {
                if (result_x) ITFAILS;
                if (result_y) ITFAILS;
            }
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Parallel equivalency passes on " << node;
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
        // >>> UNIT TESTS
        equivalence_check(ut);
        global_fails  += ut.numFails;
        global_passes += ut.numPasses;
        ut.reset();

        // sum up global errors
        nemesis::global_sum(global_fails);
        nemesis::global_sum(global_passes);
        ut.numFails  = global_fails;
        ut.numPasses = global_passes;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstParallel_Utils, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstParallel_Utils, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstParallel_Utils.cc
//---------------------------------------------------------------------------//
