//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tstLambda.cc
 * \author Thomas M. Evans
 * \date   Wed May 14 10:58:36 2014
 * \brief  Lambda testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <algorithm>

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test Helpers
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Lambda, function)
{
    std::vector<int> x = {10, 22, 31, 44, 56};

    // find first odd
    auto i = std::find_if(std::begin(x), std::end(x),
                          [](int n){ return n % 2 == 1; });

    EXPECT_EQ(31, *i);

    // find first even
    i = std::find_if(std::begin(x), std::end(x),
                     [](int n){ return n % 2 == 0; });

    EXPECT_EQ(10, *i);

    auto y = x;

    // add 1 to data
    auto f = [](int n) { return n + 1; };
    std::transform(std::begin(x), std::end(x), std::begin(y), f);

    int n = 0;
    for (const auto &m : y)
    {
        EXPECT_EQ(x[n] + 1, m);
        ++n;
    }
}

//---------------------------------------------------------------------------//
//                 end of tstLambda.cc
//---------------------------------------------------------------------------//
