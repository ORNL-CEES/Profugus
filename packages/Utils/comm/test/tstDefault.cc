//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstDefault.cc
 * \author Gregory Davidson and Stuart Slattery
 * \date   Tue Jun 21 10:52:11 2011
 * \brief  Test that we can reset the default communicator.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
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

void run_test(Parallel_Unit_Test &ut)
{
    typedef nemesis::Communicator_t         Communicator_t;

    Communicator_t new_default;
    nemesis::split(0, nemesis::nodes() - nemesis::node() - 1, new_default);

    // Get my MPI_COMM_WORLD rank
    int world_rank = nemesis::node();

    // Change the default
    nemesis::set_default(new_default);

    int new_rank = nemesis::node();

    std::cout << "World rank: " << world_rank << "   New rank: " << new_rank
              << std::endl;

    UNIT_TEST(world_rank == nemesis::nodes() - nemesis::node() - 1);

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Set default comm works on node "
          << nemesis::node() << "/" << nemesis::nodes()-1;
        ut.passes(m.str());
    }
}




//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, nemesis::release::short_version);

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        run_test(ut);
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
        std::cout << "ERROR: While testing tstDefault, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstDefault, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstDefault.cc
//---------------------------------------------------------------------------//
