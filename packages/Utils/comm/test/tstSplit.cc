//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstSplit.cc
 * \author Thomas M. Evans
 * \date   Tue Apr 13 16:00:02 2010
 * \brief  Test of comm splitting functions
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>

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

void simple_split(Parallel_Unit_Test &ut)
{
    nemesis::Communicator_t split_comm;

    // split into even/odd nodes
    if (node % 2 == 0)
        nemesis::split(0, 0, split_comm);
    else
        nemesis::split(1, 0, split_comm);

    // set the split communicator
    nemesis::set_internal_comm(split_comm);

    int new_nodes = nemesis::nodes();
    int new_node  = nemesis::node();

    if (nodes == 1)
    {
        UNIT_TEST(new_nodes == 1);
    }
    else if (nodes == 2)
    {
        UNIT_TEST(new_nodes == 1);
        UNIT_TEST(new_node == 0);
    }
    else if (nodes == 3)
    {
        if (node == 0)
        {
            UNIT_TEST(new_nodes == 2);
            UNIT_TEST(new_node == 0);
        }
        else if (node == 1)
        {
            UNIT_TEST(new_nodes == 1);
            UNIT_TEST(new_node == 0);
        }
        else if (node == 2)
        {
            UNIT_TEST(new_nodes == 2);
            UNIT_TEST(new_node == 1);
        }
    }
    else if (nodes == 4)
    {
        if (node == 0)
        {
            UNIT_TEST(new_nodes == 2);
            UNIT_TEST(new_node == 0);
        }
        else if (node == 1)
        {
            UNIT_TEST(new_nodes == 2);
            UNIT_TEST(new_node == 0);
        }
        else if (node == 2)
        {
            UNIT_TEST(new_nodes == 2);
            UNIT_TEST(new_node == 1);
        }
        else if (node == 3)
        {
            UNIT_TEST(new_nodes == 2);
            UNIT_TEST(new_node == 1);
        }
    }

    // test a reduction
    int x = 1 + 2 * node;

    nemesis::global_sum(x);

    if (nodes == 1)
    {
        UNIT_TEST(x == 1);
    }
    else if (nodes == 2)
    {
        if (node == 0)
        {
            UNIT_TEST(x == 1);
        }
        else if (node == 1)
        {
            UNIT_TEST(x == 3);
        }
    }
    else if (nodes == 3)
    {
        if (node == 0 || node == 2)
        {
            UNIT_TEST(x == 6);
        }
        else if (node == 1)
        {
            UNIT_TEST(x == 3);
        }
    }
    else if (nodes == 4)
    {
        if (node == 0 || node == 2)
        {
            UNIT_TEST(x == 6);
        }
        else if (node == 1 || node == 3)
        {
            UNIT_TEST(x == 10);
        }
    }

    // reset and destroy the split communicator
    nemesis::reset_internal_comm();
    nemesis::free_comm(split_comm);

    UNIT_TEST(nemesis::nodes() == nodes);
    UNIT_TEST(nemesis::node() == node);

    // test a reduction
    x = 1 + 2 * node;
    nemesis::global_sum(x);

    if (nodes == 1)
    {
        UNIT_TEST(x == 1);
    }
    else if (nodes == 2)
    {
        UNIT_TEST(x == 4);
    }
    else if (nodes == 3)
    {
        UNIT_TEST(x == 9);
    }
    else if (nodes == 4)
    {
        UNIT_TEST(x == 16);
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Simple splittings work on " << node;
        ut.passes(m.str());
    }
}
//---------------------------------------------------------------------------//
// the groups are:
//
//     0     1     color
// | 0 1 2 | 0 |   key
//   0 3 2   1     world rank

void four_split(Parallel_Unit_Test &ut)
{
    if (nodes != 4) return;

    nemesis::Communicator_t split_comm;
    int color = 0;
    if (node == 1) color = 1;

    // put 0, 2, and 3 in one group (in this group switch 2 and 3 in rank)
    // put 1 in another group
    if (node == 0)
        nemesis::split(color, 0, split_comm);
    else if (node == 1)
        nemesis::split(color, 0, split_comm);
    else if (node == 2)
        nemesis::split(color, 5, split_comm); // actual number doesn't matter
                                              // just the order
    else if (node == 3)
        nemesis::split(color, 3, split_comm);

    nemesis::set_internal_comm(split_comm);

    // do a broadcast from 0 on 0
    int x = 0;
    if (node == 0) x = 2;
    if (node == 1) x = 3;
    nemesis::broadcast(&x, 1, 0);

    int nnode  = nemesis::node();
    int nnodes = nemesis::nodes();

    if (node == 0)
    {
        UNIT_TEST(nnode == 0);
        UNIT_TEST(nnodes == 3);
        UNIT_TEST(x == 2);
    }
    else if (node == 1)
    {
        UNIT_TEST(nnode == 0);
        UNIT_TEST(nnodes == 1);
        UNIT_TEST(x == 3);
    }
    else if (node == 2)
    {
        UNIT_TEST(nnode == 2);
        UNIT_TEST(nnodes == 3);
        UNIT_TEST(x == 2);
    }
    else if (node == 3)
    {
        UNIT_TEST(nnode == 1);
        UNIT_TEST(nnodes == 3);
        UNIT_TEST(x == 2);
    }

    // reset and destroy the split communicator
    nemesis::reset_internal_comm();
    nemesis::free_comm(split_comm);
    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Four processor split works on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//
/*
 Reorder 0 1 2 3 with
 x =     1 2 1 1
 into    0 2 1 3
 */
void reorder_test(Parallel_Unit_Test &ut)
{
    if (nodes != 4) return;

    // make data on each node
    int n = node;
    int x = 1;
    if (node == 1)
        x = 2;

    vector<int> key(nodes, 0);
    key[node] = x;
    nemesis::global_sum(&key[0], nodes);

    int num_1  = count(key.begin(), key.end(), 1);
    int offset = 0;
    if (x == 1)
        offset = 0;
    if (x == 2)
        offset = num_1;

    int my_key = 0;
    for (int i = 0; i < node; ++i)
    {
        if (key[i] == x) my_key++;
    }

    nemesis::Communicator_t reorder_comm;

    nemesis::split(0, my_key+offset, reorder_comm);

    nemesis::set_internal_comm(reorder_comm);

    UNIT_TEST(nemesis::nodes() == 4);

    if (n == 0)
    {
        UNIT_TEST(nemesis::node() == 0);
        UNIT_TEST(x == 1);
    }
    else if (n == 1)
    {
        UNIT_TEST(nemesis::node() == 3);
        UNIT_TEST(x == 2);
    }
    else if (n == 2)
    {
        UNIT_TEST(nemesis::node() == 1);
        UNIT_TEST(x == 1);
    }
    else if (n == 3)
    {
        UNIT_TEST(nemesis::node() == 2);
        UNIT_TEST(x == 1);
    }

    // reset and destroy the split communicator
    nemesis::reset_internal_comm();
    nemesis::free_comm(reorder_comm);

    UNIT_TEST(nemesis::node() == n);

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Four processor reorder works on " << node;
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
        int gpass = 0;
        int gfail = 0;

        if (nodes < 5)
        {
            simple_split(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();

            four_split(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();

            reorder_test(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }
        else
        {
            gpass++;
        }

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstSplit, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSplit, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstSplit.cc
//---------------------------------------------------------------------------//
