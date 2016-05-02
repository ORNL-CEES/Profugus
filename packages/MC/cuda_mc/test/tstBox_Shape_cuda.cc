//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstBox_Shape_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 15:22:15 2016
 * \brief  Test for Box_Shape
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Box_Shape_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class Box_ShapeTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS

  protected:
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Box_ShapeTest, inside)
{
    std::vector<double> box_bounds = {0.0, 1.0, -1.0, 1.0, 3.0, 5.0};

    int num_vals = 4;
    std::vector<cuda_utils::Space_Vector> pts = {{1.2,   0.0, 3.5},
                                           {0.5,   0.1, 4.9},
                                           {0.75, -2.0, 4.0},
                                           {0.1,  -0.4, 5.01}};
    std::vector<int>    inside(num_vals,0);

    cuda_mc::Box_Shape_Tester::test_inside(box_bounds,pts,inside);

    std::vector<int> expected = {0, 1, 0, 0};
    EXPECT_VEC_EQ( expected, inside );
}

TEST_F(Box_ShapeTest, sample)
{
    std::vector<double> box_bounds = {0.0, 1.0, -1.0, 1.0, 3.0, 5.0};

    int num_vals = 4;
    std::vector<cuda_utils::Space_Vector> pts(num_vals);

    cuda_mc::Box_Shape_Tester::test_sample(box_bounds,pts);

    for( const auto &pt : pts )
    {
        EXPECT_TRUE( pt.x > box_bounds[0] );
        EXPECT_TRUE( pt.x < box_bounds[1] );
        EXPECT_TRUE( pt.y > box_bounds[2] );
        EXPECT_TRUE( pt.y < box_bounds[3] );
        EXPECT_TRUE( pt.z > box_bounds[4] );
        EXPECT_TRUE( pt.z < box_bounds[5] );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstBox_Shape.cc
//---------------------------------------------------------------------------//
