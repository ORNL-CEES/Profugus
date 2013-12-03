//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/test/tstRequest.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:34:36 2008
 * \brief  Non-blocking request test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "release/Release.hh"
#include "harness/Soft_Equivalence.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "../Parallel_Unit_Test.hh"

#ifdef COMM_MPI
#include "../MPI_Traits.hh"
#include <mpi.h>
#endif

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;
using nemesis::Request;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstCopyConstructor(Parallel_Unit_Test &ut)
{
    Request requestA;
    Request requestB( requestA );

    // The behavior of the copy constructor is not obvious.  If requestA has
    // not been used (inuse() returns 0) then requestA != requestB.

    if( requestA.inuse() == 0 && requestA == requestB  )
    {
        ut.passes("requestA.inuse()==0, so requestA cannot == requestB.");
    }

    if( requestA.inuse() == 0 && requestA != requestB  )
    {
        ut.passes("requestA.inuse()==0 and requestA != requestB.");
    }

    if( requestA.inuse() == 1 && requestA == requestB  )
    {
        ut.passes("requestA.inuse()=1 and requestA == requestB.");
    }

    if( requestA.inuse() == 1 && requestA != requestB  )
    {
        ut.failure("requestA.inuse()=1, so requestA must == requestB.");
    }

    if (ut.numFails == 0)
        ut.passes("Requests ok.");
}

//---------------------------------------------------------------------------//

void tstTraits(Parallel_Unit_Test &ut)
{
    using nemesis::Comm_Traits;
    {
        nemesis::HSyncSpinLock headsyncspinlock;

        if( Comm_Traits<unsigned char>::tag  != 432 ) ITFAILS;
        if( Comm_Traits<short>::tag          != 433 ) ITFAILS;
        if( Comm_Traits<unsigned short>::tag != 434 ) ITFAILS;
        if( Comm_Traits<unsigned int>::tag   != 436 ) ITFAILS;
        if( Comm_Traits<unsigned long>::tag  != 438 ) ITFAILS;
        if( Comm_Traits<long double>::tag    != 441 ) ITFAILS;

        if (ut.numFails == 0)
            ut.passes("Comm_Traits ok.");
    }
#ifdef COMM_MPI
    {
        using nemesis::MPI_Traits;

        int    i        = 2;
        short si        = 1;
        unsigned int ui = 4;
        long  li        = 8;

        float f  = 2.0;
        double d = 10.0;

        // integer type
        {
            int r = 0;
            MPI_Allreduce(&i, &r, 1, MPI_Traits<int>::element_type(), MPI_SUM,
                          MPI_COMM_WORLD);
            UNIT_TEST(r == nodes * 2);
        }

        // short integer type
        {
            short r = 0;
            MPI_Allreduce(&si, &r, 1, MPI_Traits<short>::element_type(), MPI_SUM,
                          MPI_COMM_WORLD);
            UNIT_TEST(r == nodes);
        }

        // unsigned integer type
        {
            unsigned int r = 0;
            MPI_Allreduce(&ui, &r, 1, MPI_Traits<unsigned int>::element_type(),
                          MPI_SUM, MPI_COMM_WORLD);
            UNIT_TEST(r == static_cast<unsigned int>(nodes * 4));
        }

        // long integer type
        {
            long r = 0;
            MPI_Allreduce(&li, &r, 1, MPI_Traits<long>::element_type(),
                          MPI_SUM, MPI_COMM_WORLD);
            UNIT_TEST(r == nodes * 8);
        }

        // float type
        {
            float r = 0;
            MPI_Allreduce(&f, &r, 1, MPI_Traits<float>::element_type(),
                          MPI_SUM, MPI_COMM_WORLD);
            UNIT_TEST(soft_equiv(r, static_cast<float>(nodes * 2)));
        }

        // double type
        {
            double r = 0;
            MPI_Allreduce(&d, &r, 1, MPI_Traits<double>::element_type(),
                          MPI_SUM, MPI_COMM_WORLD);
            UNIT_TEST(soft_equiv(r, static_cast<double>(nodes * 10)));
        }

        if (ut.numFails == 0)
            ut.passes("MPI_Traits ok.");
    }
#endif
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    nemesis::Parallel_Unit_Test ut(argc, argv, nemesis::release::short_version);

    node  = nemesis::node();
    nodes = nemesis::nodes();

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        tstCopyConstructor(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        tstTraits(ut);
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
        std::cout << "ERROR: While testing tstRequest, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstRequest, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstRequest.cc
//---------------------------------------------------------------------------//
