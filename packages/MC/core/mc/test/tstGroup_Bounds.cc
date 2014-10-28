//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstGroup_Bounds.cc
 * \author Thomas M. Evans
 * \date   Thu May 01 11:03:24 2014
 * \brief  Group_Bounds unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Group_Bounds.hh"

#include "gtest/utils_gtest.hh"

using profugus::Group_Bounds;
typedef Group_Bounds::Vec_Dbl Vec_Dbl;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(GroupBounds, interface)
{
    Vec_Dbl bounds(4);
    bounds[0] = 10.; bounds[1] =  5.; bounds[2] =  2.; bounds[3] =  1.;

    Group_Bounds gb(bounds);

    EXPECT_EQ(3, gb.num_groups());

    ASSERT_EQ(4, gb.group_bounds().size());
    EXPECT_DOUBLE_EQ(10., gb.group_bounds()[0]);
    EXPECT_DOUBLE_EQ( 5., gb.group_bounds()[1]);
    EXPECT_DOUBLE_EQ( 2., gb.group_bounds()[2]);
    EXPECT_DOUBLE_EQ( 1., gb.group_bounds()[3]);

    // Test finding of energies inside group bounds
    int index = -1;

    EXPECT_TRUE(gb.find(10.  , index)); EXPECT_EQ(0, index);
    EXPECT_TRUE(gb.find( 5.01, index)); EXPECT_EQ(0, index);
    EXPECT_TRUE(gb.find( 5.  , index)); EXPECT_EQ(0, index);
    EXPECT_TRUE(gb.find( 4.99, index)); EXPECT_EQ(1, index);
    EXPECT_TRUE(gb.find( 2.01, index)); EXPECT_EQ(1, index);
    EXPECT_TRUE(gb.find( 2.0 , index)); EXPECT_EQ(1, index);
    EXPECT_TRUE(gb.find( 1.99, index)); EXPECT_EQ(2, index);
    EXPECT_TRUE(gb.find( 1.5 , index)); EXPECT_EQ(2, index);
    EXPECT_TRUE(gb.find( 1.01, index)); EXPECT_EQ(2, index);
    EXPECT_TRUE(gb.find( 1.0 , index)); EXPECT_EQ(2, index);

    double upper = -1., lower = -1.;
    gb.get_energy(0, lower, upper);
    EXPECT_DOUBLE_EQ(10., upper);
    EXPECT_DOUBLE_EQ(5., lower);

    gb.get_energy(1, lower, upper);
    EXPECT_DOUBLE_EQ(5., upper);
    EXPECT_DOUBLE_EQ(2., lower);

    gb.get_energy(2, lower, upper);
    EXPECT_DOUBLE_EQ(2., upper);
    EXPECT_DOUBLE_EQ(1., lower);

    EXPECT_FALSE(gb.find( 11., index));
    EXPECT_FALSE(gb.find( 0.1, index));
}

//---------------------------------------------------------------------------//

TEST(GroupBounds, logarithmic)
{
    auto gb = Group_Bounds::build_logarithmic(
        1.0E-5, 1.0E+5, 10);

    EXPECT_EQ(10, gb->num_groups());

    double last = 1.0e6;
    for (int i = 0; i <= gb->num_groups(); ++i)
    {
        last *= 0.1;
        EXPECT_SOFTEQ(last, gb->group_bounds()[i], 1e-6);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstGroup_Bounds.cc
//---------------------------------------------------------------------------//
