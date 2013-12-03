//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/test/tstCasts.cc
 * \author Seth R Johnson
 * \date   Wed Sep 18 15:12:05 2013
 * \brief  Unit test for Casts.hh
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Casts.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Scalar_Unit_Test.hh"
#include "Release.hh"

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

using nemesis::Scalar_Unit_Test;
using nemesis::smart_cast;

//---------------------------------------------------------------------------//
// HELPER CLASSES
//---------------------------------------------------------------------------//

struct Base
{
    int base_val;

    // Base class must be polymorphic for DBC run-time checking
    virtual ~Base() {/* * */}
};

template<typename T>
struct Derived : public Base
{
    T derived_val;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_smart_cast_ref(Scalar_Unit_Test &ut)
{
    typedef Derived<float> Derived_t;
    Derived_t d;
    d.base_val = 123;
    d.derived_val = 3.5f;

    // TEST CONSTANT REFERENCE
    const Base& base_ref = d;
    UNIT_TEST(base_ref.base_val == 123);

    const Derived_t& casted_ref = smart_cast<const Derived_t&>(base_ref);
    UNIT_TEST(casted_ref.base_val == 123);
    UNIT_TEST(casted_ref.derived_val == 3.5f);

    // TEST MUTABLE REFERENCE
    Base&    base_mref = d;
    Derived_t& casted_mref = smart_cast<Derived_t&>(base_mref);
    casted_mref.derived_val = 5.0f;

    UNIT_TEST(d.derived_val == 5.0f);

    if (ut.numFails == 0)
        ut.passes("Reference smart_casts pass.");
}

void test_smart_cast_ptr(Scalar_Unit_Test &ut)
{
    typedef Derived<float> Derived_t;
    Derived_t d;
    d.base_val = 123;
    d.derived_val = 3.5f;

    // TEST CONSTANT POINTER
    const Base* base_ptr = &d;
    UNIT_TEST(base_ptr->base_val == 123);

    const Derived_t* casted_ptr = smart_cast<const Derived_t*>(base_ptr);
    UNIT_TEST(casted_ptr->base_val == 123);
    UNIT_TEST(casted_ptr->derived_val == 3.5f);

    // TEST MUTABLE POINTER
    Base*    base_mptr = &d;
    Derived_t* casted_mptr = smart_cast<Derived_t*>(base_mptr);
    casted_mptr->derived_val = 5.0f;

    UNIT_TEST(d.derived_val == 5.0f);

    // TEST NULL POINTER
    const Base* base_null_ptr = NULL;
    const Derived_t* casted_null_ptr
        = smart_cast<const Derived_t*>(base_null_ptr);

    UNIT_TEST(casted_null_ptr == NULL);

    if (ut.numFails == 0)
        ut.passes("Pointer smart_casts pass.");
}


void test_bad_smart_cast_ref(Scalar_Unit_Test &ut)
{
    typedef Derived<float>  Derived_t;

    Derived_t d;
    d.base_val = 123;
    d.derived_val = 3.5f;

    // TEST CONSTANT REFERENCE
    const Base& base_ref = d;

#ifdef REQUIRE_ON
    try
    {
        typedef Derived<double> Other_Derived_t;
        const Other_Derived_t& casted_ref
            = smart_cast<const Other_Derived_t&>(base_ref);
        // should have thrown
        ITFAILS;
        std::cout << casted_ref.derived_val << std::endl;
    }
    catch (const nemesis::assertion& a)
    {
        ut.passes("Successfully caught bad smart_cast.");
    }
    catch(...)
    {
        ut.failure("Threw the wrong type of error.");
    }

#endif

    if (ut.numFails == 0)
        ut.passes("Reference error catching pass.");
}

void test_bad_smart_cast_ptr(Scalar_Unit_Test &ut)
{
    typedef Derived<float>  Derived_t;

    Derived_t d;
    d.base_val = 123;
    d.derived_val = 3.5f;

    // TEST CONSTANT ptrERENCE
    const Base* base_ptr = &d;

#ifdef REQUIRE_ON
    try
    {
        typedef Derived<double> Other_Derived_t;
        const Other_Derived_t* casted_ptr
            = smart_cast<const Other_Derived_t*>(base_ptr);
        // should have thrown
        ITFAILS;
        std::cout << casted_ptr->derived_val << std::endl;
    }
    catch (const nemesis::assertion& a)
    {
        ut.passes("Successfully caught bad smart_cast.");
    }
    catch(...)
    {
        ut.failure("Threw the wrong type of error.");
    }
#endif

    if (ut.numFails == 0)
        ut.passes("Pointer error catching pass.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Scalar_Unit_Test ut(argc, argv, nemesis::release);

    try
    {
        test_smart_cast_ref(ut);
        test_smart_cast_ptr(ut);
        test_bad_smart_cast_ref(ut);
        test_bad_smart_cast_ptr(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstCasts, "
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstCasts, "
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                 end of tstCasts.cc
//---------------------------------------------------------------------------//
