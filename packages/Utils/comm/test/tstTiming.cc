//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstTiming.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  9 09:16:46 2008
 * \brief  Timing macros test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "release/Release.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"
#include "../Timing.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

typedef nemesis::Timing_Diagnostics D;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void do_A()
{
    TIMER(A_timer);
    TIMER_START(A_timer);

    // do a mat-vec multiply

    int S = 300000;

    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1];

    TIMER_STOP(A_timer);
    TIMER_RECORD("A_iteration", A_timer);
}

//---------------------------------------------------------------------------//

void do_B()
{
    TIMER(B_timer);
    TIMER_START(B_timer);

    // do a mat-vec multiply

    int S = 200000;

    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1];

    TIMER_STOP(B_timer);
    TIMER_RECORD("B_iteration", B_timer);
}

//---------------------------------------------------------------------------//

void do_C()
{
    TIMER(C_timer);
    TIMER_START(C_timer);

    // do a mat-vec multiply

    int S = 100000;

    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1];

    TIMER_STOP(C_timer);
    TIMER_RECORD("C_iteration", C_timer);
}


//---------------------------------------------------------------------------//

void do_D()
{
    SCOPED_TIMER("do_D");

    // do a mat-vec multiply

    int S = 100000;

    vector<double> b(S, 0.0);
    vector<double> x(S, 0.0);

    double A = 1.0;
    double B = 1.1;
    double C = 1.01;

    x[0] = 0.2;
    for (int i = 1; i < S; i++)
        x[i] = i + 1.1 + x[i-1];

    for (int i = 1; i < (S - 1); ++i)
    {
        b[i] = x[i-1] * A + x[i] * B + x[i+1] * C;
    }

    b[0]   = B * x[0] + C * x[1];
    b[S-1] = A * x[S-2] + B * x[S-1];
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void timing_active(Parallel_Unit_Test &ut)
{
#ifdef DENOVO_TIMING_ON
    cout << ">>> Testing timing macros with value " << DENOVO_TIMING << endl;
#else
    cout << ">>> Timing macros inactive" << endl;
#endif
}

//---------------------------------------------------------------------------//

void test_timing(Parallel_Unit_Test &ut)
{
    // add to some timers
    D::update_timer("A", 1.2);
    D::update_timer("B", 1.1);
    D::update_timer("B", 2.3);

    if (!soft_equiv(D::timer_value("A"), 1.2)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 3.4)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;

    D::reset_timer("B");
    D::update_timer("A", 1.3);
    if (!soft_equiv(D::timer_value("A"), 2.5)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;

    vector<string> timers = D::timer_keys();
    if (timers.size() != 3) ITFAILS;
    if (timers[0] != "A") ITFAILS;
    if (timers[1] != "B") ITFAILS;
    if (timers[2] != "C") ITFAILS;

    D::delete_timer("B");
    timers = D::timer_keys();
    if (timers.size() != 2) ITFAILS;
    if (timers[0] != "A") ITFAILS;
    if (timers[1] != "C") ITFAILS;

    // calling timer_value on B will get it back
    if (!soft_equiv(D::timer_value("A"), 2.5)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;
    timers = D::timer_keys();
    if (timers.size() != 3)   ITFAILS;
    if (D::num_timers() != 3) ITFAILS;

    // delete all timers
    D::delete_timers();
    if (D::num_timers() != 0) ITFAILS;
    timers = D::timer_keys();
    if (timers.size() != 0)   ITFAILS;

    D::update_timer("B", 12.4);
    D::update_timer("C", 1.3);
    if (!soft_equiv(D::timer_value("A"), 0.0))  ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 12.4)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 1.3))  ITFAILS;

    // reset all timers
    D::reset_timers();
    if (!soft_equiv(D::timer_value("A"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("B"), 0.0)) ITFAILS;
    if (!soft_equiv(D::timer_value("C"), 0.0)) ITFAILS;

    if (ut.numFails == 0)
    {
        ut.passes("Diagnostics timer lists ok.");
    }
}

//---------------------------------------------------------------------------//

void test_macros(Parallel_Unit_Test &ut)
{
    // delete all existing timers
    D::delete_timers();

    // make timers and do results
    TIMER(outer_timer);
    TIMER_START(outer_timer);

    do_A();
    do_B();
    do_C();
    do_D();

    TIMER_STOP(outer_timer);
    TIMER_RECORD("Outer", outer_timer);

    // if the timers are off we get no timing data
    vector<string> keys = D::timer_keys();
    if (NEMESIS_TIMING == 0)
    {
        if (keys.size() != 0)     ITFAILS;
        if (D::num_timers() != 0) ITFAILS;
    }
    else
    {
        if (keys.size() != 5) ITFAILS;
        cout << setw(15) << "Routine" << setw(15) << "Fraction" << endl;
        cout << "------------------------------" << endl;

        // get the keys and print a table
        double total = D::timer_value("Outer");
        if (total < 0.0) ITFAILS;

        if (total == 0.0)
        {
            cout << "Computer too fast for this precision timer."
                 << endl << endl;

        }
        else
        {
            cout.precision(4);
            cout.setf(ios::fixed, ios::floatfield);

            for (int i = 0, N = keys.size(); i < N; ++i)
            {
                double fraction = D::timer_value(keys[i]) / total;
                cout << setw(15) << keys[i] << setw(15) << fraction << endl;
            }

            cout << "The total time was " << total << endl;
            cout << endl;
        }
    }

    if (ut.numFails == 0)
    {
        ut.passes("Timer macros ok.");
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

        if (nodes == 1)
        {
            timing_active(ut);
            test_timing(ut);
            test_macros(ut);

            gpass = ut.numPasses;
            gfail = ut.numFails;
        }
        else
        {
            gpass = 1;
        }

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstTiming, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTiming, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstTiming.cc
//---------------------------------------------------------------------------//
