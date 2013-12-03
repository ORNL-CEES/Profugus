//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstParallel_Unit_Test.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:35:24 2008
 * \brief  Parallel_Unit_Test harness test.
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

void tstOne(nemesis::Unit_Test &unitTest)
{
    unitTest.passes("Looks like the PASSMSG macro is working.");
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
        tstOne(ut);
        ut.status();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstParallel_Unit_Test, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstParallel_Unit_Test, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstParallel_Unit_Test.cc
//---------------------------------------------------------------------------//
