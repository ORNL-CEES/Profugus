//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstReduction.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:36:47 2008
 * \brief  Global reduction tests.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iterator>

#include "harness/DBC.hh"
#include "release/Release.hh"
#include "harness/Soft_Equivalence.hh"
#include "../global.hh"
#include "../Parallel_Unit_Test.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

int node  = 0;
int nodes = 0;

using nemesis::global_sum;
using nemesis::global_prod;
using nemesis::global_min;
using nemesis::global_max;
using nemesis::sum;
using nemesis::prod;
using nemesis::min;
using nemesis::max;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void elemental_global_reduction(Parallel_Unit_Test &ut)
{
    // test ints
    int xint = nemesis::node() + 1;
    global_sum(xint);

    int int_answer = 0;
    for (int i = 0; i < nemesis::nodes(); i++)
        int_answer += i + 1;

    if (xint != int_answer) ITFAILS;

    // Test with deprecated form of global_sum
    xint = nemesis::node() + 1;
    global_sum(xint);
    if (xint != int_answer) ITFAILS;

    // test longs
    long xlong = nemesis::node() + 1000;
    global_sum(xlong);

    long long_answer = 0;
    for (int i = 0; i < nemesis::nodes(); i++)
        long_answer += i + 1000;

    if (xlong != long_answer) ITFAILS;

    // test doubles
    double xdbl = static_cast<double>(nemesis::node()) + 0.1;
    global_sum(xdbl);

    double dbl_answer = 0.0;
    for (int i = 0; i < nemesis::nodes(); i++)
        dbl_answer += static_cast<double>(i) + 0.1;

    if (!soft_equiv(xdbl, dbl_answer)) ITFAILS;

    // test product
    xlong = nemesis::node() + 1;
    global_prod(xlong);

    long_answer = 1;
    for (int i = 0; i < nemesis::nodes(); i++)
        long_answer *= (i + 1);

    if (xlong != long_answer) ITFAILS;

    // Test with deprecated form of global_prod
    xlong = nemesis::node() + 1;
    global_prod(xlong);
    if (xlong != long_answer) ITFAILS;

    // test min
    xdbl = 0.5 + nemesis::node();
    global_min(xdbl);

    if (!soft_equiv(xdbl, 0.5)) ITFAILS;

    // Test with deprecated form of global_min
    xdbl = nemesis::node() + 0.5;
    global_min(xdbl);
    if (!soft_equiv(xdbl, 0.5)) ITFAILS;

    // test max
    xdbl = 0.7 + nemesis::node();
    global_max(xdbl);

    if (!soft_equiv(xdbl, nemesis::nodes() - 0.3)) ITFAILS;

    // Test with deprecated form of global_max
    xdbl = 0.7 + nemesis::node();
    global_max(xdbl);
    if (!soft_equiv(xdbl, nemesis::nodes() - 0.3)) ITFAILS;

    if (ut.numFails == 0)
        ut.passes("Elemental global reductions ok.");
}

//---------------------------------------------------------------------------//

void array_global_reduction(Parallel_Unit_Test &ut)
{
    // make a vector of doubles
    vector<double> x(100);
    vector<double> prod(100, 1.0);
    vector<double> sum(100, 0.0);
    vector<double> min(100, 0.0);
    vector<double> max(100, 0.0);

    // fill it
    for (int i = 0; i < 100; i++)
    {
        x[i]  = nemesis::node() + 0.11;
        for (int j = 0; j < nemesis::nodes(); j++)
        {
            sum[i]  += (j + 0.11);
            prod[i] *= (j + 0.11);
        }
        min[i] = 0.11;
        max[i] = nemesis::nodes() + 0.11 - 1.0;
    }

    vector<double> c;

    {
        c = x;
        global_sum(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), sum.begin(), sum.end())) ITFAILS;

        c = x;
        global_prod(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), prod.begin(), prod.end())) ITFAILS;

        c = x;
        global_min(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), min.begin(), min.end())) ITFAILS;

        c = x;
        global_max(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), max.begin(), max.end())) ITFAILS;

    }

    // Test using deprecated forms of global_sum, global_min, global_max and
    // global_prod.

    {
        c = x;
        global_sum(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), sum.begin(), sum.end())) ITFAILS;

        c = x;
        global_prod(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), prod.begin(), prod.end())) ITFAILS;

        c = x;
        global_min(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), min.begin(), min.end())) ITFAILS;

        c = x;
        global_max(&c[0], 100);
        if (!soft_equiv(c.begin(), c.end(), max.begin(), max.end())) ITFAILS;

    }

    if (ut.numFails == 0)
        ut.passes("Array global reductions ok.");
}

//---------------------------------------------------------------------------//

void elemental_reduction(Parallel_Unit_Test &ut)
{
    // test ints
    int xint = node + 1;
    sum(xint, 0);

    int int_answer = nodes * (nodes + 1) / 2;

    if (node == 0)
    {
        UNIT_TEST(xint == int_answer);
    }
    else
    {
        UNIT_TEST(xint == node + 1);
    }

    // test doubles writing to node 2
    if (nodes == 4)
    {
        double xdbl = node + 1.0;
        sum(xdbl, 2);

        double dbl_answer = nodes * (nodes + 1.0) / 2.0;

        if (node == 2)
        {
            UNIT_TEST(soft_equiv(xdbl, dbl_answer));
        }
        else
        {
            UNIT_TEST(soft_equiv(xdbl, node + 1.0));
        }

        // test product
        xdbl       = (node + 1.5);
        dbl_answer = (1.5 * 2.5 * 3.5 * 4.5);

        prod(xdbl, 1);

        if (node == 1)
        {
            UNIT_TEST(soft_equiv(xdbl, dbl_answer));
        }
        else
        {
            UNIT_TEST(soft_equiv(xdbl, node + 1.5));
        }

        // test max
        max(xdbl, 3);

        if (node == 3 || node == 1)
        {
            UNIT_TEST(soft_equiv(xdbl, dbl_answer));
        }
        else
        {
            UNIT_TEST(soft_equiv(xdbl, node + 1.5));
        }

        // test min
        if (node == 0) xdbl = 10000.0;
        min(xdbl, 1);

        if (node == 0)
        {
            UNIT_TEST(soft_equiv(xdbl, 10000.0));
        }
        else if (node == 1)
        {
            UNIT_TEST(soft_equiv(xdbl, 3.5));
        }
        else if (node == 2)
        {
            UNIT_TEST(soft_equiv(xdbl, node + 1.5));
        }
        else
        {
            UNIT_TEST(soft_equiv(xdbl, dbl_answer));
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Elemental reductions ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//

void array_reduction(Parallel_Unit_Test &ut)
{
    if (nodes != 4) return;

    vector<long> ref(5);
    ref[0] = node + 10;
    ref[1] = node + 12;
    ref[2] = node + 13;
    ref[3] = node + 1;
    ref[4] = node + 5;

    vector<long> ans(5);
    vector<long> x(ref);

    // sum

    ans[0] = 46;
    ans[1] = 54;
    ans[2] = 58;
    ans[3] = 10;
    ans[4] = 26;

    sum(&x[0], 5, 1);

    if (node == 1)
    {
        UNIT_TEST(x[0] == ans[0]);
        UNIT_TEST(x[1] == ans[1]);
        UNIT_TEST(x[2] == ans[2]);
        UNIT_TEST(x[3] == ans[3]);
        UNIT_TEST(x[4] == ans[4]);
    }
    else
    {
        UNIT_TEST(x[0] == ref[0]);
        UNIT_TEST(x[1] == ref[1]);
        UNIT_TEST(x[2] == ref[2]);
        UNIT_TEST(x[3] == ref[3]);
        UNIT_TEST(x[4] == ref[4]);
    }

    // product
    x = ref;

    ans[0] = 17160;
    ans[1] = 32760;
    ans[2] = 43680;
    ans[3] = 24;
    ans[4] = 1680;

    prod(&x[0], 5, 3);

    if (node == 3)
    {
        UNIT_TEST(x[0] == ans[0]);
        UNIT_TEST(x[1] == ans[1]);
        UNIT_TEST(x[2] == ans[2]);
        UNIT_TEST(x[3] == ans[3]);
        UNIT_TEST(x[4] == ans[4]);
    }
    else
    {
        UNIT_TEST(x[0] == ref[0]);
        UNIT_TEST(x[1] == ref[1]);
        UNIT_TEST(x[2] == ref[2]);
        UNIT_TEST(x[3] == ref[3]);
        UNIT_TEST(x[4] == ref[4]);
    }

    // min
    x = ref;

    ans[0] = 10;
    ans[1] = 12;
    ans[2] = 13;
    ans[3] = 1;
    ans[4] = 5;

    min(&x[0], 5, 2);

    if (node == 0 || node == 2)
    {
        UNIT_TEST(x[0] == ans[0]);
        UNIT_TEST(x[1] == ans[1]);
        UNIT_TEST(x[2] == ans[2]);
        UNIT_TEST(x[3] == ans[3]);
        UNIT_TEST(x[4] == ans[4]);
    }
    else
    {
        UNIT_TEST(x[0] == ref[0]);
        UNIT_TEST(x[1] == ref[1]);
        UNIT_TEST(x[2] == ref[2]);
        UNIT_TEST(x[3] == ref[3]);
        UNIT_TEST(x[4] == ref[4]);
    }

    // max
    x = ref;

    ans[0] = 13;
    ans[1] = 15;
    ans[2] = 16;
    ans[3] = 4;
    ans[4] = 8;

    max(&x[0], 5, 1);

    if (node == 1 || node == 3)
    {
        UNIT_TEST(x[0] == ans[0]);
        UNIT_TEST(x[1] == ans[1]);
        UNIT_TEST(x[2] == ans[2]);
        UNIT_TEST(x[3] == ans[3]);
        UNIT_TEST(x[4] == ans[4]);
    }
    else
    {
        UNIT_TEST(x[0] == ref[0]);
        UNIT_TEST(x[1] == ref[1]);
        UNIT_TEST(x[2] == ref[2]);
        UNIT_TEST(x[3] == ref[3]);
        UNIT_TEST(x[4] == ref[4]);
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Array reductions ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//
/*!

  Do the following matrix-vector multiply:

  b = Ax

 */
void reduce_scatter(Parallel_Unit_Test &ut)
{
    if (nodes != 4) return;

    // reference
    vector<double>           rb(9, 0.0);
    vector<double>           rx(9, 0.0);
    vector< vector<double> > A(9, vector<double>(9, 0.0));
    {
        A[0][0] = 1.1; A[0][1] = 1.3; A[0][2] = 1.2;
        A[0][3] = 3.1; A[0][4] = 0.1; A[0][5] = 1.8;
        A[0][6] = 0.2; A[0][7] = 1.5; A[0][8] = 2.1;

        for (int i = 1; i < 9; ++i)
        {
            for (int j = 0; j < 9; ++j)
            {
                A[i][j] = A[i-1][j] + 0.6 * static_cast<double>(j+1.2);
            }
        }

        rx[0] = 1.1; rx[1] = 4.2; rx[2] = 7.8;
        rx[3] = 1.3; rx[4] = 2.5; rx[5] = 2.6;
        rx[6] = 4.1; rx[7] = 1.9; rx[8] = 4.3;

        for (int i = 0; i < 9; ++i)
        {
            for (int j = 0; j < 9; ++j)
            {
                rb[i] += A[i][j] * rx[j];
            }
        }
    }

    // steering vector
    vector<int> s(4, 0);
    s[0] = 2; s[1] = 3; s[2] = 1; s[3] = 3;

    // l2g map
    vector<int> l2g(s[node], 0);

    if (node == 0)
    {
        l2g[0] = 0;
        l2g[1] = 1;
    }

    if (node == 1)
    {
        l2g[0] = 2;
        l2g[1] = 3;
        l2g[2] = 4;
    }

    if (node == 2)
    {
        l2g[0] = 5;
    }

    if (node == 3)
    {
        l2g[0] = 6;
        l2g[1] = 7;
        l2g[2] = 8;
    }

    // local and global b
    vector<double> lb(s[node], 0.0), b(9, 0.0);

    // local x
    vector<double> x(s[node], 0.0);
    for (int i = 0; i < s[node]; ++i)
        x[i] = rx[l2g[i]];

    // on each node write the contribution into b
    for (int i = 0; i < 9; ++i)
    {
        for (int j = 0; j < s[node]; ++j)
        {
            b[i] += A[i][l2g[j]] * x[j];
        }
    }

    // reduce-scatter the final result into local b
    nemesis::sum_scatter(&b[0], &lb[0], &s[0]);

    for (int i = 0; i < s[node]; ++i)
    {
        UNIT_TEST(soft_equiv(lb[i], rb[l2g[i]]));
        UNIT_TEST(b[l2g[i]] < rb[l2g[i]]);
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Reduce-scatters ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//
/*!

  Perform all-to-all communication

 */
void all_to_all_test(Parallel_Unit_Test &ut)
{
    if (nodes != 4) return;

    int n=2;
    int total=n*nodes;
    std::vector<int> x(total);
    std::vector<int> y(total);

    // Set data to send
    x[0] = total*node;
    for (int i=0; i<total-1; ++i)
    {
        x[i+1] = x[i] + 1;
    }

    // cout << "Input on node " << node << ": ";
    // copy(x.begin(), x.end(), ostream_iterator<int>(cout, ", "));
    // cout << endl;

    // Do communication
    nemesis::all_to_all(&x[0], &y[0], n);

    // cout << "Output on node " << node << ": ";
    // copy(y.begin(), y.end(), ostream_iterator<int>(cout, ", "));
    // cout << endl;

    // Check received data
    int start = n*node;
    for (int i=0; i<nodes; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            UNIT_TEST(y[j+i*n]==start+j+i*total);
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "All-to-alls ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//
/*!

  Perform all-to-all communication with variable message size

 */
void variable_all_to_all(Parallel_Unit_Test &ut)
{
    if (nodes != 4) return;

    std::vector<double> x,y;
    std::vector<int>    sendcounts(nodes),sendoffsets(nodes);
    std::vector<int>    recvcounts(nodes),recvoffsets(nodes);

    int num_send,num_recv;
    if( node==0 )
    {
        num_send=2;
        num_recv=2;
    }
    else if( node==1 )
    {
        num_send=3;
        num_recv=0;
    }
    else if( node==2 )
    {
        num_send=0;
        num_recv=3;
    }
    else if( node==3 )
    {
        num_send=1;
        num_recv=1;
    }

    x.resize(num_send);
    y.resize(num_recv);

    if( node==0 )
    {
        x[0] = 0.1;
        x[1] = 1.2;

        sendcounts[0] = 0;
        sendcounts[1] = 0;
        sendcounts[2] = 1;
        sendcounts[3] = 1;

        sendoffsets[0] = 0;
        sendoffsets[1] = 0;
        sendoffsets[2] = 0;
        sendoffsets[3] = 1;

        recvcounts[0] = 0;
        recvcounts[1] = 1;
        recvcounts[2] = 0;
        recvcounts[3] = 1;

        recvoffsets[0] = 0;
        recvoffsets[1] = 0;
        recvoffsets[2] = 1;
        recvoffsets[3] = 1;
    }
    else if( node==1 )
    {
        x[0] = 2.3;
        x[1] = 3.4;
        x[2] = 4.5;
        y.resize(1);

        sendcounts[0] = 1;
        sendcounts[1] = 0;
        sendcounts[2] = 2;
        sendcounts[3] = 0;

        sendoffsets[0] = 0;
        sendoffsets[1] = 1;
        sendoffsets[2] = 1;
        sendoffsets[3] = 3;

        recvcounts[0] = 0;
        recvcounts[1] = 0;
        recvcounts[2] = 0;
        recvcounts[3] = 0;

        recvoffsets[0] = 0;
        recvoffsets[1] = 0;
        recvoffsets[2] = 0;
        recvoffsets[3] = 0;
    }
    else if( node==2 )
    {
        x.resize(1);

        sendcounts[0] = 0;
        sendcounts[1] = 0;
        sendcounts[2] = 0;
        sendcounts[3] = 0;

        sendoffsets[0] = 0;
        sendoffsets[1] = 0;
        sendoffsets[2] = 0;
        sendoffsets[3] = 0;

        recvcounts[0] = 1;
        recvcounts[1] = 2;
        recvcounts[2] = 0;
        recvcounts[3] = 0;

        recvoffsets[0] = 0;
        recvoffsets[1] = 1;
        recvoffsets[2] = 3;
        recvoffsets[3] = 3;
    }
    else if( node==3 )
    {
        x[0] = 5.6;

        sendcounts[0] = 1;
        sendcounts[1] = 0;
        sendcounts[2] = 0;
        sendcounts[3] = 0;

        sendoffsets[0] = 0;
        sendoffsets[1] = 1;
        sendoffsets[2] = 1;
        sendoffsets[3] = 1;

        recvcounts[0] = 1;
        recvcounts[1] = 0;
        recvcounts[2] = 0;
        recvcounts[3] = 0;

        recvoffsets[0] = 0;
        recvoffsets[1] = 1;
        recvoffsets[2] = 1;
        recvoffsets[3] = 1;
    }

    // Do communication
    nemesis::all_to_all(&x[0], &sendcounts[0], &sendoffsets[0],
                        &y[0], &recvcounts[0], &recvoffsets[0]);

    // Check received data
    if( node==0 )
    {
        UNIT_TEST(soft_equiv(y[0],2.3));
        UNIT_TEST(soft_equiv(y[1],5.6));
    }
    else if( node==2 )
    {
        UNIT_TEST(soft_equiv(y[0],0.1));
        UNIT_TEST(soft_equiv(y[1],3.4));
        UNIT_TEST(soft_equiv(y[2],4.5));
    }
    else if( node==3 )
    {
        UNIT_TEST(soft_equiv(y[0],1.2));
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Variable all-to-alls ok on " << node;
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

        elemental_global_reduction(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        array_global_reduction(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        elemental_reduction(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        array_reduction(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        reduce_scatter(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        all_to_all_test(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        variable_all_to_all(ut);
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
        std::cout << "ERROR: While testing tstReduction, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstReduction, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstReduction.cc
//---------------------------------------------------------------------------//
