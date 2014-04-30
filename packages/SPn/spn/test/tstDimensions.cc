//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstDimensions.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 23 22:38:34 2012
 * \brief  SPN Dimensions test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Dimensions.hh"

using profugus::Dimensions;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Dimensions, SP1)
{
    Dimensions dim(1);

    EXPECT_EQ(1, dim.num_equations());
    EXPECT_EQ(2, dim.num_moments());
    EXPECT_EQ(4, Dimensions::max_num_equations());
}

//---------------------------------------------------------------------------//

TEST(Dimensions, SP3)
{
    Dimensions dim(3);

    EXPECT_EQ(2, dim.num_equations());
    EXPECT_EQ(4, dim.num_moments());
}

//---------------------------------------------------------------------------//

TEST(Dimensions, SP5)
{
    Dimensions dim(5);

    EXPECT_EQ(3, dim.num_equations());
    EXPECT_EQ(6, dim.num_moments());
}

//---------------------------------------------------------------------------//

TEST(Dimensions, SP7)
{
    Dimensions dim(7);

    EXPECT_EQ(4, dim.num_equations());
    EXPECT_EQ(8, dim.num_moments());
}

//---------------------------------------------------------------------------//
//                        end of tstDimensions.cc
//---------------------------------------------------------------------------//
