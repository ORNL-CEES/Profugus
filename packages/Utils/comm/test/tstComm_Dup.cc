//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstComm_Dup.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:34:59 2008
 * \brief  Communicator duplication test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "release/Release.hh"
#include "../SpinLock.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_mpi_comm_dup(Parallel_Unit_Test &ut)
{
    // we only run this particular test when mpi is on
#ifdef COMM_MPI

    Require (nemesis::nodes() == 4);

    int node  = nemesis::node();
    int snode = 0;

    // split up nodes (two communicators) 0 -> 0, 2 -> 1 and
    // 1 -> 0, 3 -> 1
    MPI_Comm new_comm;

    if (node == 1)
    {
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &new_comm);
        MPI_Comm_rank(new_comm, &snode);
    }
    else if (node == 3)
    {
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &new_comm);
        MPI_Comm_rank(new_comm, &snode);
    }
    else
    {
        MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &new_comm);
        MPI_Comm_rank(new_comm, &snode);
    }

    // we haven't set the communicator yet so we should still have 4 nodes
    if (nemesis::nodes() != 4) ITFAILS;

    // now dup the communicator on each processor
    nemesis::inherit(new_comm);

    // each processor should see two nodes
    if (nemesis::nodes() != 2) ITFAILS;

    // test data send/receive
    int data = 0;

    // do some tests on each processor
    if (node == 0)
    {
        if (nemesis::node() != 0) ITFAILS;

        // set data to 10 and send it out
        data = 10;
        nemesis::send(&data, 1, 1, 100);
    }
    else if (node == 1)
    {
        if (nemesis::node() != 0) ITFAILS;

        // set data to 20 and send it out
        data = 20;
        nemesis::send(&data, 1, 1, 100);
    }
    else if (node == 2)
    {
        if (nemesis::node() != 1) ITFAILS;

        if (data != 0) ITFAILS;
        nemesis::receive(&data, 1, 0, 100);
        if (data != 10) ITFAILS;
    }
    else if (node == 3)
    {
        if (nemesis::node() != 1) ITFAILS;

        if (data != 0) ITFAILS;
        nemesis::receive(&data, 1, 0, 100);
        if (data != 20) ITFAILS;
    }

    // now free the communicator on each processor
    nemesis::free_inherited_comm();

    // the free should have set back to COMM_WORLD
    if (nemesis::nodes() != 4) ITFAILS;

    {
        nemesis::HTSyncSpinLock slock;
        for (int i = 0; i < 10000; i++)
        {
            continue;
        }

        if (ut.numFails == 0)
        {
            ostringstream m;
            m << "Communicator duplicated successfully on " << nemesis::node();
            ut.passes(m.str());
        }
    }

    nemesis::global_barrier();

    MPI_Comm_free(&new_comm);

#endif
}

//---------------------------------------------------------------------------//

void test_comm_dup(Parallel_Unit_Test &ut)
{
    // we only run this test scalar
#ifdef COMM_SCALAR

    int node = nemesis::node();

    // now dup the communicator on each processor
    nemesis::inherit(node);

    // check the number of nodes
    if (nemesis::nodes() != 1) ITFAILS;

    nemesis::free_inherited_comm();

    // check the number of nodes
    if (nemesis::nodes() != 1) ITFAILS;

    if (ut.numFails == 0)
        ut.passes("Scalar Comm duplication/free works ok.");

#endif

    // check duping/freeing MPI_COMM_WORLD
#ifdef COMM_MPI

    int nodes = nemesis::nodes();

    MPI_Comm comm_world = MPI_COMM_WORLD;
    nemesis::inherit(comm_world);

    if (nemesis::nodes() != nodes) ITFAILS;

    // try a global sum to check
    int x = 10;
    nemesis::global_sum(x);
    if (x != 10 * nodes) ITFAILS;

    nemesis::free_inherited_comm();

    // we should be back to comm world
    if (nemesis::nodes() != nodes) ITFAILS;

    // try a global sum to check
    int y = 20;
    nemesis::global_sum(y);
    if (y != 20 * nodes) ITFAILS;

    if (ut.numFails == 0)
        ut.passes("MPI_COMM_WORLD Comm duplication/free works ok.");

#endif
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

        if (nemesis::nodes() == 4)
        {
            test_mpi_comm_dup(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }

        test_comm_dup(ut);
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
        std::cout << "ERROR: While testing tstComm_Dup, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstComm_Dup, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstComm_Dup.cc
//---------------------------------------------------------------------------//
