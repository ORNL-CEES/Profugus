//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Array_Tester.hh
 * \author Tom Evans
 * \date   Wed Jan 04 23:29:15 2017
 * \brief  RTK_Array_Tester class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_test_RTK_Array_Tester_hh
#define MC_cuda_rtk_test_RTK_Array_Tester_hh

#include <memory>

#include "Utils/gtest/Gtest_Functions.hh"
#include "MC/geometry/RTK_Cell.hh"
#include "MC/geometry/RTK_Array.hh"

//---------------------------------------------------------------------------//
// TESTERS
//---------------------------------------------------------------------------//

class Base : public ::testing::Test
{
  protected:
    using Lattice = profugus::RTK_Array<profugus::RTK_Cell>;
    using Core    = profugus::RTK_Array<Lattice>;
    using SP_Core = std::shared_ptr<Core>;

  protected:
    virtual ~Base() = default;

  protected:
    SP_Core core;
};

//---------------------------------------------------------------------------//

class SimpleLattice : public Base
{
  protected:

    void SetUp()
    {
        using SP_Object = Lattice::SP_Object;
        using Object_t  = Lattice::Object_t;

        // make a 3x2x1 lattice with 10 pin-types (only 3 defined)
        lattice = std::make_shared<Lattice>(3, 2, 1, 10);
        Lattice& lat = *lattice;

        // pins 0, 3 and 5
        SP_Object pin0(std::make_shared<Object_t>(0, 0.62, 10, 1.26, 14.28));
        SP_Object pin3(std::make_shared<Object_t>(3, 0.45, 10, 1.26, 14.28));
        SP_Object pin5(std::make_shared<Object_t>(5, 0.54, 10, 1.26, 14.28));

        CHECK(2 == pin0->num_cells());
        CHECK(2 == pin3->num_cells());
        CHECK(2 == pin5->num_cells());

        // assign pins to the lattice
        lat.id(0, 0, 0) = 3;
        lat.id(1, 0, 0) = 3;
        lat.id(2, 0, 0) = 3;
        lat.id(0, 1, 0) = 5;
        lat.id(1, 1, 0) = 0;
        lat.id(2, 1, 0) = 5;

        // assign pins to ids
        lat.assign_object(pin0, 0);
        lat.assign_object(pin3, 3);
        lat.assign_object(pin5, 5);

        lat.complete(0.0, 0.0, 0.0);

        run_test();
    }

    void run_test();

  protected:
    std::shared_ptr<Lattice> lattice;
};

//---------------------------------------------------------------------------//

class SimpleCore : public Base
{
  protected:

    void SetUp()
    {
        core = std::make_shared<Core>(3, 3, 4, 1);

        // 1 pin
        Lattice::SP_Object pin(std::make_shared<Lattice::Object_t>(
                                   1, 0.5, 5, 1.5, 4.0));
        EXPECT_EQ(2, pin->num_cells());

        // 1 lattice
        Core::SP_Object lat(std::make_shared<Lattice>(2, 2, 1, 1));
        lat->assign_object(pin,  0);
        lat->id(0, 0, 0) = 0;
        lat->id(1, 0, 0) = 0;
        lat->id(0, 1, 0) = 0;
        lat->id(1, 1, 0) = 0;
        lat->complete(0.0, 0.0, 0.0);

        // assign lattices to core
        core->assign_object(lat, 0);
        for (int k = 0; k < 4; ++k)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int i = 0; i < 3; ++i)
                {
                    core->id(i, j, k) = 0;
                }
            }
        }

        core->complete(1.1, 0.0, -5.0);

        run_test();
    }

    void run_test();
};

//---------------------------------------------------------------------------//

class RegCore : public Base
{
  protected:

    void SetUp()
    {
        // 3x3 core with 3 objects (2 lattice types, one reflector) + water on
        // top
        core = std::make_shared<Core>(3, 3, 2, 3);
        {
            // 2 pin types

            // local id = 0
            Lattice::SP_Object black = std::make_shared<Lattice::Object_t>(
                1, 0.54, 5, 1.26, 14.28);
            EXPECT_EQ(2, black->num_cells());

            // local id = 0
            Lattice::SP_Object purple = std::make_shared<Lattice::Object_t>(
                2, 0.54, 5, 1.26, 14.28);
            EXPECT_EQ(2, purple->num_cells());

            // local id = 0
            Lattice::SP_Object water = std::make_shared<Lattice::Object_t>(
                5, 2.52, 14.28);
            EXPECT_EQ(1, water->num_cells());

            // 3 lattices
            Core::SP_Object lat1 = std::make_shared<Lattice>(2, 2, 1, 1);
            Core::SP_Object lat2 = std::make_shared<Lattice>(2, 2, 1, 1);
            Core::SP_Object lat0 = std::make_shared<Lattice>(1, 1, 1, 1);

            // assign pins to ids
            lat1->assign_object(black,  0);
            lat2->assign_object(purple, 0);
            lat0->assign_object(water, 0);

            // assign pins in the lattices
            lat1->id(0, 0, 0) = 0;
            lat1->id(1, 0, 0) = 0;
            lat1->id(0, 1, 0) = 0;
            lat1->id(1, 1, 0) = 0;

            lat2->id(0, 0, 0) = 0;
            lat2->id(1, 0, 0) = 0;
            lat2->id(0, 1, 0) = 0;
            lat2->id(1, 1, 0) = 0;

            lat0->id(0, 0, 0) = 0;

            lat1->complete(0.0, 0.0, 0.0);
            lat2->complete(0.0, 0.0, 0.0);
            lat0->complete(0.0, 0.0, 0.0);

            // assign lattices to core
            core->assign_object(lat1, 1);
            core->assign_object(lat2, 2);
            core->assign_object(lat0, 0);

            // assign ids
            core->id(0, 0, 0) = 1;
            core->id(1, 0, 0) = 2;
            core->id(2, 0, 0) = 0;
            core->id(0, 1, 0) = 2;
            core->id(1, 1, 0) = 1;
            core->id(2, 1, 0) = 0;
            core->id(0, 2, 0) = 0;
            core->id(1, 2, 0) = 0;
            core->id(2, 2, 0) = 0;
            core->id(0, 0, 1) = 0;
            core->id(1, 0, 1) = 0;
            core->id(2, 0, 1) = 0;
            core->id(0, 1, 1) = 0;
            core->id(1, 1, 1) = 0;
            core->id(2, 1, 1) = 0;
            core->id(0, 2, 1) = 0;
            core->id(1, 2, 1) = 0;
            core->id(2, 2, 1) = 0;
        }
    }

    void run_test();
    void reflect_test();
};

#endif // MC_cuda_rtk_test_RTK_Array_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Array_Tester.hh
//---------------------------------------------------------------------------//
