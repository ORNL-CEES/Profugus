//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstDemangle_Type.cc
 * \author Seth R Johnson
 * \date   Wed Oct 16 08:35:37 2013
 * \brief  Demangle_Type test.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Demangle_Type.hh"

#include <iostream>

#include "harness/Scalar_Unit_Test.hh"
#include "Release.hh"

using namespace std;
using namespace nemesis;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

struct A
{
    virtual ~A()
    {
    }
};

struct B : public A
{
    int unused;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_demangle(Scalar_Unit_Test &ut)
{
    // >>> Simple types
    std::string int_type = Demangle_Type::demangle(typeid(int).name());
    std::string flt_type = Demangle_Type::demangle(typeid(float).name());

    cout << "Demangled int: '"   << int_type << "'" << endl;
    cout << "Demangled float: '" << flt_type << "'" << endl;

    UNIT_TEST(int_type != flt_type);

#ifdef __GNUG__
    UNIT_TEST(int_type == "int");
#endif

    // >>> Class types
    std::string sut_type = Demangle_Type::demangled_type<Scalar_Unit_Test>();
    std::string str_type = Demangle_Type::demangled_type<std::string>();

    cout << "Demangled Scalar_Unit_Test: '" << sut_type << "'" << endl;
    cout << "Demangled std::string: '"      << str_type << "'" << endl;

#ifdef __GNUG__
    UNIT_TEST(sut_type == "nemesis::Scalar_Unit_Test");
#endif

    UNIT_TEST(sut_type != str_type);

    // Make sure that the same call gives us the same answer
    std::string sut_type2 = Demangle_Type::demangled_type<Scalar_Unit_Test>();

    UNIT_TEST(sut_type == sut_type2);

    // >>> Dynamic typename
    std::string a_type = Demangle_Type::demangle(typeid(A).name());
    std::string b_type = Demangle_Type::demangle(typeid(B).name());
    A* temp = new B;
    std::string temp_type = Demangle_Type::demangled_type(*temp);
    delete temp;

    cout << "Demangled A: '"    << a_type    << "'" << endl;
    cout << "Demangled B: '"    << b_type    << "'" << endl;
    cout << "Demangled temp: '" << temp_type << "'" << endl;

    UNIT_TEST(temp_type == b_type);

    if (ut.numFails == 0)
        ut.passes("Demangling works.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Scalar_Unit_Test ut(argc, argv, release);

    try
    {
        // >>> UNIT TESTS
        test_demangle(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstDemangle_Type, "
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstDemangle_Type, "
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstDemangle_Type.cc
//---------------------------------------------------------------------------//
