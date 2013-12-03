//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstComm_FC.cc
 * \author Thoams M. Evans
 * \date   Thu Sep  3 23:36:40 2009
 * \brief  FORTRAN interface test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <comm/config.h>

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
// F90 TEST FUNCTION DECLARATIONS
//---------------------------------------------------------------------------//

#define BUILD_COMM FC_FUNC_(build_comm, BUILD_COMM)
#define TEST_RANK FC_FUNC_(test_rank, TEST_RANK)
#define TEST_BROADCAST FC_FUNC_(test_broadcast, TEST_BROADCAST)
#define TEST_REDUCTION FC_FUNC_(test_reduction, TEST_REDUCTION)
#define TEST_BLOCKING FC_FUNC_(test_blocking, TEST_BLOCKING)

extern "C"
{
    void BUILD_COMM();
    void TEST_RANK(unsigned int *, unsigned int *, int *, int *);
    void TEST_BROADCAST(unsigned int *, unsigned int *);
    void TEST_REDUCTION(unsigned int *, unsigned int *);
    void TEST_BLOCKING(unsigned int *, unsigned int *);
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void basic_tests(Parallel_Unit_Test &ut)
{
    unsigned int pass = 0, fail = 0;

    // call FORTRAN simple test
    TEST_RANK(&pass, &fail, &node, &nodes);
    UNIT_TEST(fail == 0);

    if (fail == 0)
    {
        ostringstream m;
        m << "FORTRAN rank ok on " << node;
        ut.passes(m.str());
    }
    ut.numPasses += pass;
    ut.numFails  += fail;
    pass = 0; fail = 0;

    // broadcast test
    TEST_BROADCAST(&pass, &fail);
    UNIT_TEST(fail == 0);

    if (fail == 0)
    {
        ostringstream m;
        m << "FORTRAN broadcasts ok on " << node;
        ut.passes(m.str());
    }
    ut.numPasses += pass;
    ut.numFails  += fail;
    pass = 0; fail = 0;

    // reduction test
    TEST_REDUCTION(&pass, &fail);
    UNIT_TEST(fail == 0);

    if (fail == 0)
    {
        ostringstream m;
        m << "FORTRAN global reductions ok on " << node;
        ut.passes(m.str());
    }
    ut.numPasses += pass;
    ut.numFails  += fail;
    pass = 0; fail = 0;
}

//---------------------------------------------------------------------------//

void blocking_comm(Parallel_Unit_Test &ut)
{
    unsigned int pass = 0, fail = 0;

    // call FORTRAN simple test
    TEST_BLOCKING(&pass, &fail);
    UNIT_TEST(fail == 0);

    if (fail == 0)
    {
        ostringstream m;
        m << "FORTRAN blocking send/receive ok on " << node;
        ut.passes(m.str());
    }
    ut.numPasses += pass;
    ut.numFails  += fail;
    pass = 0; fail = 0;
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

#if MPI_IMP != MPI_MPT

        // build the FORTRAN comm module
        BUILD_COMM();

        // tests
        basic_tests(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        if (nodes == 2)
        {
            blocking_comm(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }

#else

        ostringstream m;
        m << "FORTRAN Comm interface requires MPI 2.0 standard; mpt < 2.0! "
          << "FORTRAN interface unavailable.";
        ut.passes(m.str());
        gpass++;

#endif

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstComm_FC, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstComm_FC, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstComm_FC.cc
//---------------------------------------------------------------------------//
