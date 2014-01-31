//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstStatic_Map.cc
 * \author Seth R Johnson
 * \date   Tue Dec 17 21:04:10 2013
 * \brief  Static_Map test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Static_Map.hh"

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <map>

using std::size_t;
using profugus::Static_Map;

//---------------------------------------------------------------------------//
// TEST FIXTURE
//---------------------------------------------------------------------------//
class StaticMap : public ::testing::Test
{
  protected:
    typedef Static_Map<size_t, double> Map_t;
    typedef Map_t::value_type          value_type; // key/val pair
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(StaticMap, zero_elements)
{
    // make a hash table that needs to hold 0 elements
    Map_t ht;

    EXPECT_FALSE(ht.completed());
    ht.complete();

    EXPECT_EQ(251, ht.bucket_count());
    EXPECT_EQ(0, ht.size());
    EXPECT_TRUE(ht.empty());
    EXPECT_TRUE(ht.completed());

    EXPECT_TRUE(!ht.exists(1002));
    EXPECT_TRUE(!ht.exists(1));
    EXPECT_TRUE(!ht.exists(3));
    EXPECT_TRUE(!ht.exists(25251));
    EXPECT_TRUE(!ht.exists(343));
    EXPECT_TRUE(!ht.exists(0));

    EXPECT_EQ(0, ht.count(343));

    EXPECT_EQ(0, ht.size());
}

//---------------------------------------------------------------------------//
TEST_F(StaticMap, duplicate_key)
{
    // make a hash table that needs to hold 0 elements
    Map_t ht;

    ht.insert(value_type(3, 2.));
    ht.insert(value_type(3, 4.));

    EXPECT_THROW(ht.complete(), profugus::assertion);
}

//---------------------------------------------------------------------------//
TEST_F(StaticMap, 10k_elements)
{
    Map_t ht;

    for (std::size_t i = 0; i < 10000; ++i)
    {
        ht.insert(Map_t::value_type(i, i * 2));
    }

    ht.complete();

    EXPECT_EQ(16381, ht.bucket_count());
    EXPECT_EQ(10000, ht.size());

    EXPECT_TRUE(ht.exists(1001));
    EXPECT_TRUE(ht.exists(0));
    EXPECT_TRUE(ht.exists(2));
    EXPECT_TRUE(ht.exists(9165));
    EXPECT_TRUE(ht.exists(342));

    EXPECT_EQ(1, ht.count(342));
    EXPECT_EQ(0, ht.count(10000));

    EXPECT_TRUE(!ht.exists(10001));
    EXPECT_TRUE(!ht.exists(25251));
    EXPECT_TRUE(ht.find(25251) == ht.end());

    // check them
    EXPECT_EQ(2002, ht[1001]);
    EXPECT_EQ(0,    ht[0]);
    EXPECT_EQ(4,    ht[2]);
    EXPECT_EQ(684,  ht[342]);

    EXPECT_EQ(1001, ht.find(1001)->first);
    EXPECT_EQ(2002, ht.find(1001)->second);

    // check them
    EXPECT_EQ(2002, ht.at(1001));

    // change value
    ht[2003] = 2.4;

    EXPECT_EQ(2.4, ht[2003]);

#ifdef REQUIRE_ON
    // try accessing a non-inserted key
    EXPECT_THROW(ht[10001] = 2.5, profugus::assertion);

    // try inserting over an already inserted pair
    EXPECT_THROW(ht.insert(value_type(1001, 1.7)), profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
TEST_F(StaticMap, collisions)
{
    Map_t ht;

    EXPECT_EQ(0, ht.size());

    // add values (with at least one collision)
    ht.insert(value_type(21, 5.1));
    ht.insert(value_type(1025, 5.6)); // collision at bucket 21

    ht.complete();

    EXPECT_EQ(251, ht.bucket_count());
    EXPECT_EQ(2, ht.size());
    EXPECT_EQ(2, ht.bucket_size(21));
    EXPECT_EQ(21, ht.hash_function().hash(1025));

    ht.clear();

    ht.insert(value_type(21, 5.1));
    ht.insert(value_type(1025, 5.6)); // collision at bucket 21
    ht.insert(value_type(2782, 5.7)); // another collision at bucket 21

    ht.complete();

    EXPECT_EQ(3, ht.size());
    EXPECT_EQ(3, ht.bucket_size(21));
    EXPECT_EQ(3, ht.max_items_per_bucket());
}

//---------------------------------------------------------------------------//
TEST_F(StaticMap, constness)
{
    Map_t ht;

    EXPECT_EQ(0, ht.size());

    // insert some key,value objects
    ht.insert(value_type(1001,  1.1));
    ht.insert(value_type(0,     1.2));
    ht.insert(value_type(2,     1.3));
    ht.insert(value_type(25250, 1.4));
    ht.insert(value_type(342,   1.5));

    ht.complete();

    // make a const-reference
    {
        const Static_Map<size_t, double> &h = ht;

        EXPECT_TRUE(h.exists(1001));
        EXPECT_TRUE(h.exists(0));
        EXPECT_TRUE(h.exists(2));
        EXPECT_TRUE(h.exists(25250));
        EXPECT_TRUE(h.exists(342));

        // check them
        EXPECT_EQ(1.1, h[1001]);

        // check them
        EXPECT_EQ(1.1, h.at(1001));
        EXPECT_EQ(1.2, h.at(0));
        EXPECT_EQ(1.3, h.at(2));
        EXPECT_EQ(1.4, h.at(25250));
        EXPECT_EQ(1.5, h.at(342));

        EXPECT_EQ(5, h.size());
        EXPECT_FALSE(h.empty());
    }
}

//---------------------------------------------------------------------------//
class StaticMapSigned : public ::testing::Test
{
  protected:
    typedef Static_Map<int, double> Map_t;
    typedef Map_t::value_type       value_type; // key/val pair
};

//---------------------------------------------------------------------------//
TEST_F(StaticMapSigned, basic_checks)
{
    Map_t ht;

    EXPECT_EQ(0, ht.size());

    // insert some key,value objects
    ht.insert(value_type(1001, 1.1));
    ht.insert(value_type(0,    2.1));
    ht.insert(value_type(-1,   3.1));

    ht.complete();
    EXPECT_EQ(3, ht.size());

    EXPECT_TRUE(ht.exists(1001));
    EXPECT_TRUE(ht.exists(0));
    EXPECT_TRUE(ht.exists(-1));

    // check them
    EXPECT_EQ(1.1, ht.at(1001));
    EXPECT_EQ(2.1, ht.at(0));
    EXPECT_EQ(3.1, ht.at(-1));
}

//---------------------------------------------------------------------------//
TEST_F(StaticMap, map_conversion)
{
    std::map<int, double> building;

    building[4]  = 100;
    building[-3] = 98;
    building[3]  = 50;

    for (std::size_t i = 100; i < 1000; ++i)
    {
        building[i * 2] = i;
    }

    // Build hash table from map
    Map_t ht(building.begin(), building.end());
    ht.complete();

    ht[3] = 1.234;
    EXPECT_EQ(1.234, ht[3]);

    EXPECT_EQ(100., ht[4] );
    EXPECT_EQ(98.,  ht[-3]);

    for (std::size_t i = 100; i < 1000; ++i)
    {
        EXPECT_EQ(i, ht[i * 2]);
    }

    EXPECT_EQ(3 + 900, ht.size());

    // Build a map from the hash table
    std::map<int, double> checking(ht.begin(), ht.end());

    EXPECT_EQ(3 + 900, checking.size());

    EXPECT_EQ(100., checking[4] );
    EXPECT_EQ(98.,  checking[-3]);
    EXPECT_EQ(1.234, checking[3]);

    for (std::size_t i = 100; i < 1000; ++i)
    {
        EXPECT_EQ(i, checking[i * 2]);
    }
}

//---------------------------------------------------------------------------//
//                        end of tstStatic_Map.cc
//---------------------------------------------------------------------------//
