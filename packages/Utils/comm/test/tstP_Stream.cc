//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstP_Stream.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 26 10:22:20 2013
 * \brief  Test of Parallel Streaming classes.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
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
#include "../SpinLock.hh"
#include "../P_Stream.hh"

using nemesis::pout;
using nemesis::pnout;
using nemesis::pcout;
using nemesis::pncout;
using nemesis::endl;
using nemesis::fixed;
using nemesis::scientific;
using nemesis::setw;
using nemesis::setfill;
using nemesis::setprecision;

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);
#define EXPECT_EQ(a, b) if (a != b) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void Native_pcout(Parallel_Unit_Test &ut)
{
    EXPECT_EQ(0, pcout.master());
    pcout << endl;
    pcout << "Parallel cout on " << nemesis::nodes()
          << " nodes from node " << nemesis::node() << endl;

    double x = 12.45122435;
    std::cout.precision(2);
    pcout << "x = " << setw(12) << fixed << x << endl;
    pcout << "x = " << setw(12) << scientific << x << endl;

    pcout << "x = " << setprecision(5) << setw(15) << fixed << x << endl;
    pcout << "x = " << setprecision(5) << setw(15) << scientific << x << endl;

    std::cout.precision(5);
    pcout << "x = " << setfill('.') << setw(15) << fixed << x << endl;
    pcout << "x = " << setfill('.') << setw(15) << scientific << x << endl;

    pcout << endl;

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Native pcout test ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

void Native_pncout(Parallel_Unit_Test &ut)
{
    EXPECT_EQ(nemesis::node(), pncout.master());
    pcout << endl;
    nemesis::global_barrier();

    {
        nemesis::HTSyncSpinLock l;
        pncout << "Parallel-node cout on" << nemesis::nodes()
               << " nodes from node " << nemesis::node() << endl;
    }

    nemesis::global_barrier();
    pcout << endl;

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Native pncout test ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

void Local(Parallel_Unit_Test &ut)
{
    nemesis::P_Out p(1);

    EXPECT_EQ(1, p.master());
    p << endl;
    p << "Local parallel out on  " << nemesis::nodes()
      << " nodes from node " << nemesis::node() << endl;
    p << "\n";

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Local test ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

void Toggle_pout(Parallel_Unit_Test &ut)
{
    pcout << endl;
    pcout << "No message below indicates that pout is toggled off." << endl;
    pout << "Nemesis toggle is on; output on node " << nemesis::node()
         << " is active" << endl;
    pcout << endl;

#ifdef NEMESIS_POUT
    EXPECT_EQ(0, pout.master());
#else
    EXPECT_EQ(-1, pout.master());
#endif

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Toggle pout test ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

void Toggle_pnout(Parallel_Unit_Test &ut)
{
    pcout << endl;
    pcout << "No message below indicates that pnout is toggled off." << endl;
    nemesis::global_barrier();

    {
        nemesis::HTSyncSpinLock l;
        pnout << "Nemesis toggle is on; output on node " << nemesis::node()
              << " is active" << endl;
    }

    nemesis::global_barrier();
    pcout << endl;

#ifdef NEMESIS_POUT
    EXPECT_EQ(nemesis::node(), pnout.master());
#else
    EXPECT_EQ(-1, pnout.master());
#endif

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Toggle pnout test ok on " << node;
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

        Native_pcout(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        Native_pncout(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        Local(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        Toggle_pout(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        Toggle_pnout(ut);
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
        std::cout << "ERROR: While testing tstP_Stream, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstP_Stream, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstP_Stream.cc
//---------------------------------------------------------------------------//
