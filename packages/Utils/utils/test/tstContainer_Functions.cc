//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/test/tstContainer_Functions.cc
 * \author Seth R Johnson
 * \date   Sat Sep 15 16:12:23 2012
 * \brief  Tests the functions in the Container_Functions file.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <algorithm>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <functional>

#include "utils/Constants.hh"
#include "utils/Definitions.hh"
#include "utils/Container_Functions.hh"

//---------------------------------------------------------------------------//
TEST(RemoveSubset, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;
    typedef def::Vec_String Vec_String;

    // Make dbl_vec1 (superset of dbl_vec2)
    Vec_Dbl dbl_vec1;
    dbl_vec1.push_back(1.0);
    dbl_vec1.push_back(2.0);
    dbl_vec1.push_back(3.0);
    dbl_vec1.push_back(4.0);
    dbl_vec1.push_back(5.0);

    // Make dbl_vec2 (subset of dbl_vec1)
    Vec_Dbl dbl_vec2;
    dbl_vec2.push_back(1.0);
    dbl_vec2.push_back(3.00000000001);
    dbl_vec2.push_back(5.0);

    // Examine result
    profugus::remove_subset(dbl_vec1, dbl_vec2,
                            profugus::SoftIsEqual<double>(1.0e-6));
    EXPECT_TRUE(dbl_vec1.size() == 2);
    EXPECT_TRUE( soft_equiv(dbl_vec1[0], 2.0) );
    EXPECT_TRUE( soft_equiv(dbl_vec1[1], 4.0) );

    // Make str_vec1 (superset of str_vec2)
    Vec_String str_vec1;
    str_vec1.push_back("aaa");
    str_vec1.push_back("bbb");
    str_vec1.push_back("ccc");
    str_vec1.push_back("ddd");
    str_vec1.push_back("eee");

    // Make str_vec2 (subset of str_vec1)
    Vec_String str_vec2;
    str_vec2.push_back("aaa");
    str_vec2.push_back("ccc");
    str_vec2.push_back("eee");

    // Examine result
    profugus::remove_subset(str_vec1, str_vec2);
    EXPECT_TRUE(str_vec1.size() == 2);
    EXPECT_TRUE(str_vec1[0] == "bbb");
    EXPECT_TRUE(str_vec1[1] == "ddd");
}

//---------------------------------------------------------------------------//
TEST(ContainerMerge, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;

    // Make vec1 {0.0, 1.1, 2.2, 3.3, 4.4, 5.5}
    Vec_Dbl vec1(6, 0.0);
    for (unsigned int i = 0; i < vec1.size(); ++i)
    {
        vec1[i] = i * 1.1;
    }

    // Make vec2 {1.3, 3.3, 6.5}
    Vec_Dbl vec2(3);
    vec2[0] = 1.3; vec2[1] = 3.3; vec2[2] = 6.5;

    // Make the reference union result {0.0, 1.1, 1.3, 2.2, 3.3, 4.4, 5.5, 6.5}
    Vec_Dbl ref1(8, 0.0);
    ref1[1] = 1.1; ref1[2] = 1.3; ref1[3] = 2.2; ref1[4] = 3.3;
    ref1[5] = 4.4; ref1[6] = 5.5; ref1[7] = 6.5;

    // Merge vec1 and vec2 together
    profugus::container_merge(vec1, vec2, 0.001);
    EXPECT_TRUE(vec1.size() == 8);
    for (unsigned int i = 0; i < vec1.size(); ++i)
    {
        EXPECT_SOFTEQ(vec1[i], ref1[i], 1.0e-6);
    }

    // Try merging vec1 with vec3 (an empty vector)
    Vec_Dbl vec3;
    profugus::container_merge(vec3, vec1, 0.001);
    // Result should be identical to vec1
    EXPECT_TRUE(vec3.size() == vec1.size());
    for (unsigned int i = 0; i < vec1.size(); ++i)
    {
        EXPECT_SOFTEQ(vec3[i], vec1[i], 1.0e-6);
    }

    // Try merging a vector that has some elements very close to vec1
    Vec_Dbl vec4(4, 0.0);
    vec4[0] = 0.5;  vec4[1] = 1.100001;  vec4[2] = 3.5;  vec4[3] = 4.400001;
    profugus::container_merge(vec4, vec1, 0.0001);
    EXPECT_TRUE(vec4.size() == vec1.size() + 2);
    Vec_Dbl ref2(9, 0.0);
    ref2[0] = 0.0;  ref2[1] = 0.5;  ref2[2] = 1.1;  ref2[3] = 1.3;
    ref2[4] = 2.2;  ref2[5] = 3.3;  ref2[6] = 3.5;  ref2[7] = 4.4;
    ref2[8] = 5.5;
    for(unsigned int i = 0; i < 9; ++i)
    {
        EXPECT_SOFTEQ(vec4[i], ref2[i], 1.0e-6);
    }
}

//---------------------------------------------------------------------------//
// Test associative copy
TEST(Copy, all)
{
    typedef std::map<int, double>      Map_Int_Dbl;
    typedef std::multimap<int, double> MultiMap_Int_Dbl;
    typedef def::Vec_Int               Vec_Int;
    typedef def::Vec_Dbl               Vec_Dbl;
    typedef std::list<int>             List_Int;
    typedef std::list<double>          List_Dbl;

    // Create the map and multimap
    Map_Int_Dbl map;
    map.insert( std::make_pair(0, 4.0) );
    map.insert( std::make_pair(2, 5.0) );
    map.insert( std::make_pair(4, 6.0) );
    map.insert( std::make_pair(6, 7.0) );

    MultiMap_Int_Dbl multimap;
    multimap.insert( std::make_pair(0, 4.0) );
    multimap.insert( std::make_pair(0, 5.0) );
    multimap.insert( std::make_pair(1, 6.0) );
    multimap.insert( std::make_pair(1, 7.0) );

    // First test -- copy map into vectors
    {
        // Copy key
        Vec_Int vec_int(4);
        profugus::copy(map.begin(), map.end(), vec_int.begin(),
                       profugus::CopyType::KEY);
        EXPECT_EQ(vec_int[0], 0);
        EXPECT_EQ(vec_int[1], 2);
        EXPECT_EQ(vec_int[2], 4);
        EXPECT_EQ(vec_int[3], 6);

        // Copy values
        Vec_Dbl vec_dbl(4);
        profugus::copy(map.begin(), map.end(), vec_dbl.begin(),
                       profugus::CopyType::VALUE);
        EXPECT_SOFTEQ(vec_dbl[0], 4.0, 1.0e-12);
        EXPECT_SOFTEQ(vec_dbl[1], 5.0, 1.0e-12);
        EXPECT_SOFTEQ(vec_dbl[2], 6.0, 1.0e-12);
        EXPECT_SOFTEQ(vec_dbl[3], 7.0, 1.0e-12);
    }

    // Second test -- copy map into list
    {
        // Copy key
        List_Int list_int(4);
        profugus::copy(map.begin(), map.end(), list_int.begin(),
                       profugus::CopyType::KEY);
        List_Int::const_iterator list_int_iter = list_int.begin();
        EXPECT_EQ(*list_int_iter, 0);
        ++list_int_iter;
        EXPECT_EQ(*list_int_iter, 2);
        ++list_int_iter;
        EXPECT_EQ(*list_int_iter, 4);
        ++list_int_iter;
        EXPECT_EQ(*list_int_iter, 6);

        // Copy values
        List_Dbl list_dbl(4);
        profugus::copy(map.begin(), map.end(), list_dbl.begin(),
                       profugus::CopyType::VALUE);
        List_Dbl::const_iterator list_dbl_iter = list_dbl.begin();
        EXPECT_SOFTEQ(*list_dbl_iter, 4.0, 1.0e-12);
        ++list_dbl_iter;
        EXPECT_SOFTEQ(*list_dbl_iter, 5.0, 1.0e-12);
        ++list_dbl_iter;
        EXPECT_SOFTEQ(*list_dbl_iter, 6.0, 1.0e-12);
        ++list_dbl_iter;
        EXPECT_SOFTEQ(*list_dbl_iter, 7.0, 1.0e-12);
    }

    // Third test -- copy multimap into vectors
    {
        // Copy key
        Vec_Int vec_int(4);
        profugus::copy(multimap.begin(), multimap.end(), vec_int.begin(),
                       profugus::CopyType::KEY);
        EXPECT_EQ(vec_int[0], 0);
        EXPECT_EQ(vec_int[1], 0);
        EXPECT_EQ(vec_int[2], 1);
        EXPECT_EQ(vec_int[3], 1);

        // Copy values
        Vec_Dbl vec_dbl(4);
        profugus::copy(multimap.begin(), multimap.end(), vec_dbl.begin(),
                       profugus::CopyType::VALUE);
        EXPECT_SOFTEQ(vec_dbl[0], 4.0, 1.0e-12);
        EXPECT_SOFTEQ(vec_dbl[1], 5.0, 1.0e-12);
        EXPECT_SOFTEQ(vec_dbl[2], 6.0, 1.0e-12);
        EXPECT_SOFTEQ(vec_dbl[3], 7.0, 1.0e-12);
    }

    // Fourth test -- copy multimap into list
    {
        // Copy key
        List_Int list_int(4);
        profugus::copy(multimap.begin(), multimap.end(), list_int.begin(),
                       profugus::CopyType::KEY);
        List_Int::const_iterator list_int_iter = list_int.begin();
        EXPECT_EQ(*list_int_iter, 0);
        ++list_int_iter;
        EXPECT_EQ(*list_int_iter, 0);
        ++list_int_iter;
        EXPECT_EQ(*list_int_iter, 1);
        ++list_int_iter;
        EXPECT_EQ(*list_int_iter, 1);

        // Copy values
        List_Dbl list_dbl(4);
        profugus::copy(multimap.begin(), multimap.end(), list_dbl.begin(),
                       profugus::CopyType::VALUE);
        List_Dbl::const_iterator list_dbl_iter = list_dbl.begin();
        EXPECT_SOFTEQ(*list_dbl_iter, 4.0, 1.0e-12);
        ++list_dbl_iter;
        EXPECT_SOFTEQ(*list_dbl_iter, 5.0, 1.0e-12);
        ++list_dbl_iter;
        EXPECT_SOFTEQ(*list_dbl_iter, 6.0, 1.0e-12);
        ++list_dbl_iter;
        EXPECT_SOFTEQ(*list_dbl_iter, 7.0, 1.0e-12);
    }
}

//---------------------------------------------------------------------------//
// Fill map
TEST(FillMap, all)
{
    typedef std::deque<int> Deque_Int;
    typedef def::Vec_Dbl    Vec_Dbl;

    // Create a deque of keys
    Deque_Int keys;
    keys.push_back(0);
    keys.push_back(2);
    keys.push_back(4);
    keys.push_back(6);

    // Create a vector of values
    Vec_Dbl values;
    values.push_back(10.0);
    values.push_back(100.0);
    values.push_back(1000.0);
    values.push_back(10000.0);

    // Create a map
    std::map<int, double> test_map;
    profugus::fill_map(keys.begin(), keys.end(), values.begin(), test_map);

    EXPECT_EQ(test_map.size(), 4);
    // Test keys
    EXPECT_EQ(test_map.count(0), 1);
    EXPECT_EQ(test_map.count(2), 1);
    EXPECT_EQ(test_map.count(4), 1);
    EXPECT_EQ(test_map.count(6), 1);
    // Test values
    EXPECT_SOFTEQ(test_map.find(0)->second, 10.0, 1.0e-12);
    EXPECT_SOFTEQ(test_map.find(2)->second, 100.0, 1.0e-12);
    EXPECT_SOFTEQ(test_map.find(4)->second, 1000.0, 1.0e-12);
    EXPECT_SOFTEQ(test_map.find(6)->second, 10000.0, 1.0e-12);
}

//---------------------------------------------------------------------------//
// Fill multimap
TEST(FillMultiMap, all)
{
    typedef std::deque<int> Deque_Int;
    typedef def::Vec_Dbl    Vec_Dbl;

    // Create a deque of keys
    Deque_Int keys;
    keys.push_back(0);
    keys.push_back(0);
    keys.push_back(1);
    keys.push_back(1);

    // Create a vector of values
    Vec_Dbl values;
    values.push_back(10.0);
    values.push_back(100.0);
    values.push_back(1000.0);
    values.push_back(10000.0);

    // Create a map
    std::multimap<int, double> test_map;
    profugus::fill_multimap(keys.begin(), keys.end(), values.begin(), test_map);

    EXPECT_EQ(test_map.size(), 4);
    // Test keys
    EXPECT_EQ(test_map.count(0), 2);
    EXPECT_EQ(test_map.count(1), 2);
    // Test values
    EXPECT_SOFTEQ(test_map.find(0)->second, 10.0, 1.0e-12);
    EXPECT_SOFTEQ(test_map.find(1)->second, 1000.0, 1.0e-12);
}

//---------------------------------------------------------------------------//
// The the split map function
TEST(SplitMap, all)
{
    typedef std::map<int, double>      Map_Int_Dbl;
    typedef def::Vec_Int               Vec_Int;
    typedef def::Vec_Dbl               Vec_Dbl;

    // Create the map
    Map_Int_Dbl map_cont;
    map_cont.insert( std::make_pair(0, 4.0) );
    map_cont.insert( std::make_pair(2, 5.0) );
    map_cont.insert( std::make_pair(4, 6.0) );
    map_cont.insert( std::make_pair(6, 7.0) );

    // Create vectors to hold keys and values
    Vec_Int keys(4);
    Vec_Dbl values(4);
    profugus::split_map(map_cont.begin(), map_cont.end(),
                      keys.begin(), values.begin());

    // Check keys
    EXPECT_EQ(keys.size(), 4);
    EXPECT_EQ(keys[0], 0);
    EXPECT_EQ(keys[1], 2);
    EXPECT_EQ(keys[2], 4);
    EXPECT_EQ(keys[3], 6);

    // Check values
    EXPECT_EQ(values.size(), 4);
    EXPECT_SOFTEQ(values[0], 4.0, 1.0e-12);
    EXPECT_SOFTEQ(values[1], 5.0, 1.0e-12);
    EXPECT_SOFTEQ(values[2], 6.0, 1.0e-12);
    EXPECT_SOFTEQ(values[3], 7.0, 1.0e-12);
}

//---------------------------------------------------------------------------//
// Test trim
TEST(Trim, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;

    // Create a vector holding several values
    Vec_Dbl vec;
    vec.push_back(0.3);
    vec.push_back(0.5);
    vec.push_back(0.7);
    vec.push_back(1.0);
    vec.push_back(1.3);
    vec.push_back(1.5);

    // Trim the values from 0.7 to 1.3
    Vec_Dbl new_vec_1 = vec;
    profugus::trim(new_vec_1, 0.7, 1.3);

    EXPECT_TRUE(new_vec_1.size() == 3);
    EXPECT_SOFTEQ(new_vec_1[0], 0.7, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_1[1], 1.0, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_1[2], 1.3, 1.0e-6);

    // Trim the values from 0.0 to 1.3
    Vec_Dbl new_vec_2 = vec;
    profugus::trim(new_vec_2, 0.0, 1.3);
    EXPECT_TRUE(new_vec_2.size() == 5);
    EXPECT_SOFTEQ(new_vec_2[0], 0.3, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_2[1], 0.5, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_2[2], 0.7, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_2[3], 1.0, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_2[4], 1.3, 1.0e-6);

    // Trim the values from 0.9 to 2.0
    Vec_Dbl new_vec_3 = vec;
    profugus::trim(new_vec_3, 0.9, 2.0);
    EXPECT_TRUE(new_vec_3.size() == 3);
    EXPECT_SOFTEQ(new_vec_3[0], 1.0, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_3[1], 1.3, 1.0e-6);
    EXPECT_SOFTEQ(new_vec_3[2], 1.5, 1.0e-6);
}

//---------------------------------------------------------------------------//
// Test make_unique
TEST(MakeUnique, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;
    typedef def::Vec_Int    Vec_Int;
    typedef def::Vec_String Vec_Str;

    // >>> VEC_DBL TEST
    // Create a non-sorted array of doubles with some repeats
    double dbl_test1[] = {1.0, 3.0, 2.5, 2.0, 2.5, 3.0, 4.0, 8.0, 1.0};
    // Create the reference array of unique doubles
    double dbl_ref_test1[] = {1.0, 2.0, 2.5, 3.0, 4.0, 8.0};
    // Wrap with vectors
    Vec_Dbl dbl_test1_v(&dbl_test1[0], &dbl_test1[0]+9);
    Vec_Dbl dbl_test1_ref_v(&dbl_ref_test1[0], &dbl_ref_test1[0]+6);

    // Run test
    profugus::make_unique(dbl_test1_v);
    EXPECT_EQ(dbl_test1_v.size(), dbl_test1_ref_v.size());
    for(unsigned int i = 0; i < dbl_test1_ref_v.size(); ++i)
    {
        EXPECT_SOFTEQ(dbl_test1_ref_v[i], dbl_test1_v[i], 1.0e-6);
    }

    // >>> VEC_INT TEST
    // Create a non-sorted array of ints
    int int_test2[] = {1, 3, 2, 4, 2, 3, 8, 4, 5};
    // Create the reference array of unique ints
    int int_ref_test2[] = {1, 2, 3, 4, 5, 8};
    // Wrap with ints
    Vec_Int int_test2_v(&int_test2[0], &int_test2[0]+9);
    Vec_Int int_test2_ref_v(&int_ref_test2[0], &int_ref_test2[0]+6);

    // Run test
    profugus::make_unique(int_test2_v);
    EXPECT_EQ(int_test2_v.size(), int_test2_ref_v.size());
    for(unsigned int i = 0; i < int_test2_ref_v.size(); ++i)
    {
        EXPECT_SOFTEQ(int_test2_ref_v[i], int_test2_v[i], 1.0e-6);
    }

    // >>> Vec_Str TEST
    // Make a vector of strings with repeats
    Vec_Str str_test;
    str_test.push_back("aaa");
    str_test.push_back("bbb");
    str_test.push_back("ccc");
    str_test.push_back("aaa");
    str_test.push_back("ddd");
    str_test.push_back("bbb");
    str_test.push_back("ccc");

    // Make the reference vector
    Vec_Str str_ref;
    str_ref.push_back("aaa");
    str_ref.push_back("bbb");
    str_ref.push_back("ccc");
    str_ref.push_back("ddd");

    // Get and test the result
    profugus::make_unique(str_test);
    EXPECT_EQ(str_test.size(), str_ref.size() );
    for(unsigned int i = 0; i < str_test.size(); ++i)
    {
        EXPECT_EQ(str_test[i], str_ref[i]);
    }
}

//---------------------------------------------------------------------------//
// Test the arithmetic mean
TEST(CalcMean, all)
{
    typedef def::Vec_Dbl   Vec_Dbl;
    typedef std::list<int> List_Int;

    // Test the mean of a vector of doubles
    double test1[] = {2.0, 4.0, 6.0, 8.0};
    Vec_Dbl test_1v(&test1[0], &test1[0]+4);
    EXPECT_SOFTEQ(profugus::mean(test_1v.begin(), test_1v.end()), 5.0, 1.0e-6);

    // Test the mean of a list of ints
    int test2[] = {2, 4, 6, 8};
    List_Int test_2l(&test2[0], &test2[0]+4);
    EXPECT_EQ(profugus::mean(test_2l.begin(), test_2l.end()), 5);
}

//---------------------------------------------------------------------------//
// Test the arithmetic median
TEST(CalcMedian, all)
{
    typedef def::Vec_Dbl        Vec_Dbl;
    typedef std::list<int>      List_Int;

    // Test the mean of a vector of doubles
    double test1[] = {2.0, 4.0, 6.0, 8.0};
    Vec_Dbl test_1v(&test1[0], &test1[0]+4);
    EXPECT_SOFTEQ(profugus::median(test_1v.begin(), test_1v.end()), 5.0,
                  1.0e-6);

    // Test the mean of a list of ints
    int test2[] = {9, 2, 4, 6, 8};
    List_Int test_2l(&test2[0], &test2[0]+5);
    EXPECT_EQ(profugus::median(test_2l.begin(), test_2l.end()), 6);
}

//---------------------------------------------------------------------------//

class TruncateAndNormTest : public ::testing::Test
{
  protected:
    typedef std::vector<double> Vec_Dbl;

  protected:
    void SetUp()
    {
        v.push_back(0.1);
        v.push_back(0.4);
        v.push_back(0.25);
        v.push_back(0.2);
        v.push_back(0.05);
    }
  protected:
    Vec_Dbl v;
};

TEST_F(TruncateAndNormTest, beg)
{
    double result = profugus::truncate_and_norm(v.begin(), v.begin()+2, v);
    EXPECT_SOFT_EQ(result, 0.5);

    double expected[] = {.2, .8};
    EXPECT_VEC_SOFT_EQ(expected, v);
}

TEST_F(TruncateAndNormTest, end)
{
    double result = profugus::truncate_and_norm(v.begin()+2, v.begin()+5, v);
    EXPECT_SOFT_EQ(result, 0.5);

    double expected[] = {.5, .4, .1};
    EXPECT_VEC_SOFT_EQ(expected, v);
}

TEST_F(TruncateAndNormTest, degenerate)
{
    double result = profugus::truncate_and_norm(v.begin(), v.begin()+1, v);
    EXPECT_SOFT_EQ(result, 0.1);

    double expected[] = {1.};
    EXPECT_VEC_SOFT_EQ(expected, v);
}

TEST_F(TruncateAndNormTest, odd_normalization)
{
    for (Vec_Dbl::iterator it = v.begin(); it != v.end(); ++it)
        *it *= 32;

    double result = profugus::truncate_and_norm(v.begin(), v.begin()+2, v);
    EXPECT_SOFT_EQ(result, 32 * 0.5);

    double expected[] = {.2, .8};
    EXPECT_VEC_SOFT_EQ(expected, v);
}

//---------------------------------------------------------------------------//
//                        end of tstContainer_Functions.cc
//---------------------------------------------------------------------------//
