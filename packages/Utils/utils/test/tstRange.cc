//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/test/tstRange.cc
 * \author Seth R Johnson
 * \date   Sun Sep 20 10:03:34 2015
 * \brief  Range class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Range.hh"

#include <cstdint>
#include "Utils/gtest/utils_gtest.hh"

using profugus::range;
using profugus::count;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(RangeTest, ints)
{
    typedef std::vector<int> Vec_Int;

    Vec_Int vals;
    for (auto i : range(0, 4))
    {
        ASSERT_EQ(sizeof(int), sizeof(i));
        vals.push_back(i);
    }
    EXPECT_VEC_EQ((Vec_Int{0,1,2,3}), vals);
}

TEST(RangeTest, chars)
{
    for (auto i : range('A', 'Z'))
    {
        cout << i;
    }
    cout << endl;
}

TEST(RangeTest, uint_range)
{
    typedef std::vector<unsigned int> Vec_UInt;

    Vec_UInt vals;

    for (auto u : range(20u, 25u).step(2u))
    {
        vals.push_back(u);
    }

    EXPECT_VEC_EQ((Vec_UInt{20,22,24}), vals);
}

TEST(RangeTest, large)
{
    using large_int = std::uint_least64_t;

    // Note: you can't pass both 0 (int) and large_int(10) , because the range
    // can't figure out T
    for (auto i : range<large_int>(0, 10))
    {
        ASSERT_EQ(sizeof(large_int), sizeof(i));
    }
}

TEST(RangeTest, just_end)
{
    typedef std::vector<int> Vec_Int;
    Vec_Int vals;

    for (auto i : range(4))
    {
        vals.push_back(i);
        ASSERT_LT(i, 10);
    }
    EXPECT_VEC_EQ((Vec_Int{0,1,2,3}), vals);
}

TEST(RangeTest, count_default)
{
    typedef std::vector<int> Vec_Int;
    Vec_Int vals;

    auto counter = count<int>().begin();
    vals.push_back(*counter++);
    vals.push_back(*counter++);
    vals.push_back(*counter++);

    EXPECT_VEC_EQ((Vec_Int{0,1,2}), vals);
}

TEST(RangeTest, count)
{
    typedef std::vector<int> Vec_Int;
    Vec_Int vals;

    for (auto i : count(10).step(15))
    {
        if (i > 90)
            break;
        vals.push_back(i);
    }
    EXPECT_VEC_EQ((Vec_Int{10,25,40,55,70,85}), vals);
}

TEST(RangeTest, vec_fill)
{
    typedef std::vector<int> Vec_Int;
    auto r = profugus::range(0, 5);
    Vec_Int vals(r.begin(), r.end());

    EXPECT_VEC_EQ((Vec_Int{0,1,2,3,4}), vals);

    // Re-assign
    vals.assign(r.begin(), r.end());
    EXPECT_VEC_EQ((Vec_Int{0,1,2,3,4}), vals);
}

//---------------------------------------------------------------------------//
// end of Utils/utils/test/tstRange.cc
//---------------------------------------------------------------------------//
