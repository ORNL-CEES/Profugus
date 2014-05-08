//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cxx11/tstLoops.cc
 * \author Thomas M. Evans
 * \date   Wed May 07 22:49:15 2014
 * \brief  C++-11 Loop testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Loops, begin_end)
{
    int x[] = {1, 2, 3, 4};

    int n = 1;
    for (auto i = std::begin(x); i != std::end(x); ++i)
    {
        EXPECT_EQ(n, *i);
        ++n;
    }
    EXPECT_EQ(5, n);
}

//---------------------------------------------------------------------------//

TEST(Loops, const_range)
{
    int x[] = {1, 2, 3, 4};
    std::vector<int> y(3);
    y[0] = 10; y[1] = 11; y[2] = 12;

    int n = 1;
    for (auto i : x)
    {
        EXPECT_EQ(n, i);
        ++n;
    }
    EXPECT_EQ(5, n);

    n = 10;
    for (auto i : y)
    {
        EXPECT_EQ(n, i);
        ++n;
    }
    EXPECT_EQ(13, n);
}

//---------------------------------------------------------------------------//

TEST(Loops, mutable_range)
{
    int x[] = {1, 2, 3, 4};
    std::vector<int> y(3);
    y[0] = 10; y[1] = 11; y[2] = 12;

    for (auto &i : x)
    {
        i += 10;
    }

    for (auto i : y)
    {
        i += 10;
    }

    EXPECT_EQ(11, x[0]);
    EXPECT_EQ(12, x[1]);
    EXPECT_EQ(13, x[2]);
    EXPECT_EQ(14, x[3]);

    EXPECT_EQ(10, y[0]);
    EXPECT_EQ(11, y[1]);
    EXPECT_EQ(12, y[2]);

    for (auto &i : y)
    {
        i += 10;
    }

    EXPECT_EQ(20, y[0]);
    EXPECT_EQ(21, y[1]);
    EXPECT_EQ(22, y[2]);

    for (int i = 0; i < 3; ++i)
    {
        y[i] += 1;
    }

    EXPECT_EQ(21, y[0]);
    EXPECT_EQ(22, y[1]);
    EXPECT_EQ(23, y[2]);
}

//---------------------------------------------------------------------------//
//                 end of tstLoops.cc
//---------------------------------------------------------------------------//
