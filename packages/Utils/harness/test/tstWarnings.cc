//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstWarnings.cc
 * \author Thomas M. Evans
 * \date   Sun Feb 26 21:50:58 2012
 * \brief  Warnings unit test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "../DBC.hh"
#include "../Scalar_Unit_Test.hh"
#include "../Soft_Equivalence.hh"
#include "Release.hh"
#include "../Warnings.hh"

using namespace std;
using namespace nemesis;
using namespace nemesis;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void warnings_test(Scalar_Unit_Test &ut)
{
    // add some warnings
    ADD_WARNING("A");
    ADD_WARNING("B");
    ADD_WARNING("C" << "D" << " is a warning on " << 0);

    double d = 0.452452542;
    if (true)
        ADD_WARNING("this double is " << scientific << setw(16) << d);
    else
        ITFAILS;

    (void)sizeof(d); // hide unused variable warning

    // manually add warning (this is always added, regardless whether
    // NEMESIS_WARNINGS is true or not)
    NEMESIS_WARNINGS.add("Always give this warning.");

#ifdef NEMESIS_WARNINGS_ENABLED

#   ifdef NEMESIS_WARNINGS_IMMEDIATE
    // warnings should have been printed immediately
    UNIT_TEST(NEMESIS_WARNINGS.num_warnings() == 1);
#   else

    UNIT_TEST(NEMESIS_WARNINGS.num_warnings() == 5);

    //while (!NEMESIS_WARNINGS.empty())
    //{
    //    cout << NEMESIS_WARNINGS.pop() << endl;
    //}

    UNIT_TEST(!NEMESIS_WARNINGS.empty());

    UNIT_TEST(NEMESIS_WARNINGS.pop() == "A");
    UNIT_TEST(NEMESIS_WARNINGS.pop() == "B");
    UNIT_TEST(NEMESIS_WARNINGS.pop() == "CD is a warning on 0");
    NEMESIS_WARNINGS.pop();
    UNIT_TEST(NEMESIS_WARNINGS.pop()
            == "Always give this warning.");

    UNIT_TEST(NEMESIS_WARNINGS.empty());
#   endif

    if (ut.numFails == 0)
        ut.passes("Warnings correctly reported.");

#else

    UNIT_TEST(NEMESIS_WARNINGS.num_warnings() == 1);

    while (!NEMESIS_WARNINGS.empty())
    {
        cout << NEMESIS_WARNINGS.pop() << endl;
    }

    if (ut.numFails == 0)
        ut.passes("Automatic warnings correctly turned off.");

#endif
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Scalar_Unit_Test ut(argc, argv, nemesis::release);

    try
    {
        // >>> UNIT TESTS
        warnings_test(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstWarnings, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstWarnings, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstWarnings.cc
//---------------------------------------------------------------------------//
