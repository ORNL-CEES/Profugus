//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstTypes.cc
 * \author Thomas M. Evans
 * \date   Fri Sep  4 10:37:09 2009
 * \brief  Test of C/C++/FORTRAN types.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <typeinfo>

#include <harness/config.h>

#include "../DBC.hh"
#include "../Scalar_Unit_Test.hh"
#include "../Soft_Equivalence.hh"
#include "Release.hh"

using namespace std;
using namespace nemesis;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void sizes_test(Scalar_Unit_Test &ut)
{
#ifdef NEMESIS_INT4_SET
    {
        NEMESIS_INT4 x = 0;
        UNIT_TEST(sizeof(x) == 4);

        ostringstream m;
        m << "4-byte (32-bit) integers defined correctly as "
          << typeid(NEMESIS_INT4).name();
        ut.passes(m.str());
    }
#else
    {
        ut.fails("4-byte (32-bit) integers undefined.");
    }
#endif

#ifdef NEMESIS_INT8_SET
    {
        NEMESIS_INT8 x = 0;
        UNIT_TEST(sizeof(x) == 8);

        ostringstream m;
        m << "8-byte (64-bit) integers defined correctly as "
          << typeid(NEMESIS_INT8).name();
        ut.passes(m.str());
    }
#endif

#ifdef NEMESIS_REAL4_SET
    {
        NEMESIS_REAL4 x = 0.0;
        UNIT_TEST(sizeof(x) == 4);

        ostringstream m;
        m << "4-byte (32-bit) reals defined correctly as "
          << typeid(NEMESIS_REAL4).name();
        ut.passes(m.str());
    }
#else
    {
        ut.fails("4-byte (32-bit) reals undefined.");
    }
#endif

#ifdef NEMESIS_REAL8_SET
    {
        NEMESIS_REAL8 x = 0.0;
        UNIT_TEST(sizeof(x) == 8);

        ostringstream m;
        m << "8-byte (64-bit) reals defined correctly as "
          << typeid(NEMESIS_REAL8).name();
        ut.passes(m.str());
    }
#else
    {
        ut.fails("8-byte (64-bit) reals undefined.");
    }
#endif
}

//---------------------------------------------------------------------------//
// NOTE: Fortran testing is disabled.
#if 0

#define TEST_ARRAYS FC_FUNC_(test_arrays, TEST_ARRAYS)
extern "C"
{
    void TEST_ARRAYS(NEMESIS_INT4 *, NEMESIS_REAL4 *, NEMESIS_REAL8 *, int *);
}

void CPP_F_array_test(Scalar_Unit_Test &ut)
{
    vector<NEMESIS_INT4>  i(4, 0);
    vector<NEMESIS_REAL4> f(3, 0.0);
    vector<NEMESIS_REAL8> d(3, 0.0);

    i[0] = 2; i[1] = 4; i[2] = 6; i[3] = 8;
    f[0] = 2.1; f[1] = 4.2; f[2] = 6.3;
    d[0] = 5.2; d[1] = 6.4; d[2] = 7.6;

    int fails = 0;
    TEST_ARRAYS(&i[0], &f[0], &d[0], &fails);
    UNIT_TEST(!fails);

    if (ut.numFails == 0)
        ut.passes("C++ types successfully passed to FORTRAN.");
}

#else
void CPP_F_array_test(Scalar_Unit_Test &ut) { /* * */ }
#endif

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Scalar_Unit_Test ut(argc, argv, nemesis::release);

    try
    {
        // >>> UNIT TESTS
        sizes_test(ut);
        CPP_F_array_test(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstTypes, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTypes, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstTypes.cc
//---------------------------------------------------------------------------//
