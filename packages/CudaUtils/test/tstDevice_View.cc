//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/tstDevice_View.cc
 * \author Tom Evans
 * \date   Fri Nov 18 11:54:56 2016
 * \brief  Tests for class Device_View.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "gtest/utils_gtest.hh"

#include "Device_View_Tester.hh"

//---------------------------------------------------------------------------//
// FIXTURE
//---------------------------------------------------------------------------//

class DeviceViewTest : public ::testing::Test
{
  protected:
    using Monkey = host::Monkey;
    using Cage   = host::Cage;
    using Zoo    = host::Zoo;
    using SP_Zoo = std::shared_ptr<Zoo>;

  protected:

    void SetUp()
    {
        // Make some monkey cages
        Cage::Monkeys group1 = {std::make_shared<Monkey>(1, 3.1),
                                std::make_shared<Monkey>(2, 2.8),
                                std::make_shared<Monkey>(3, 2.9)};
        Cage::Monkeys group2 = {std::make_shared<Monkey>(4, 2.1),
                                std::make_shared<Monkey>(5, 2.4)};

        // Make the cages
        auto cage1 = std::make_shared<Cage>(group1);
        auto cage2 = std::make_shared<Cage>(group2);
        Zoo::Cages cages = {cage1, cage2};

        // Make the zoo
        zoo = std::make_shared<Zoo>(cages);
    }

  protected:

    SP_Zoo zoo;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(DeviceViewTest, one_level)
{
    test_one_level(*zoo);
}

//---------------------------------------------------------------------------//

TEST_F(DeviceViewTest, two_level)
{
    test_two_level(zoo);
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/tstDevice_View.cc
//---------------------------------------------------------------------------//
