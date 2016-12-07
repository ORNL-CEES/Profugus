//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Cell_Tester.hh
 * \author Tom Evans
 * \date   Tue Nov 29 17:08:48 2016
 * \brief  RTK_Cell_Tester class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_test_RTK_Cell_Tester_hh
#define MC_cuda_rtk_test_RTK_Cell_Tester_hh

#include <memory>
#include <vector>

#include "Utils/gtest/Gtest_Functions.hh"
#include "MC/geometry/RTK_Cell.hh"

//---------------------------------------------------------------------------//
// TESTERS
//---------------------------------------------------------------------------//

class Base : public ::testing::Test
{
  protected:
    using RTK_Cell = profugus::RTK_Cell;
    using SP_Cell  = std::shared_ptr<RTK_Cell>;
    using Vec_Cell = std::vector<SP_Cell>;

  protected:
    virtual ~Base() = default;

  protected:
    Vec_Cell pins;
};

//---------------------------------------------------------------------------//

class Single_Shell : public Base
{
  protected:

    void SetUp()
    {
        SP_Cell pin1 = std::make_shared<RTK_Cell>(1, 0.54, 10, 1.26, 14.28);
        SP_Cell pin2 = std::make_shared<RTK_Cell>(1, 0.45, 2, 1.2, 14.28);
        pins         = {pin1, pin2};
    }

    void construct();
    void track();
};

//---------------------------------------------------------------------------//

class Multi_Shell : public Base
{
  protected:

    void SetUp()
    {
        // make pin with clad
        std::vector<int>    ids = {1, 2};
        std::vector<double> rad = {0.49, 0.54};
        SP_Cell pin1 = std::make_shared<RTK_Cell>(ids, rad, 3, 1.26, 14.28);

        ids = {1, 2, 5};
        rad = {0.27, 0.486, 0.54};
        SP_Cell pin2 = std::make_shared<RTK_Cell>(ids, rad, 10, 1.26, 14.28, 4);

        pins = {pin1, pin2};
    }

    void construct();
    void track();

    void multiseg_construct();
    void multiseg_track();
};

//---------------------------------------------------------------------------//

class Empty : public Base
{
  protected:

    void SetUp()
    {
        SP_Cell pin1 = std::make_shared<RTK_Cell>(11, 1.26, 14.28);
        SP_Cell pin2 = std::make_shared<RTK_Cell>(11, 0.25, 1.26, 14.28);
        pins         = {pin1, pin2};
    }

    void square();
    void rectangle();
};

//---------------------------------------------------------------------------//

#endif // MC_cuda_rtk_test_RTK_Cell_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Cell_Tester.hh
//---------------------------------------------------------------------------//
