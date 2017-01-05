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

#endif // MC_cuda_rtk_test_RTK_Array_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Array_Tester.hh
//---------------------------------------------------------------------------//
