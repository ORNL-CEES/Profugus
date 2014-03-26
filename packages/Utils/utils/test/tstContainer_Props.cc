//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstContainer_Props.cc
 * \author Seth R Johnson
 * \date   Sat Sep 15 16:12:23 2012
 * \brief  Tests the functions in the Container_Functions file.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
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
#include "utils/Container_Props.hh"

//---------------------------------------------------------------------------//
TEST(IsSorted, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;
    typedef def::Vec_Int    Vec_Int;
    typedef def::Vec_String Vec_String;

    // Test a vector of sorted doubles
    Vec_Dbl vector_1;    // {-2.0, -1.5, 0.2, 0.8, 1.3}
    vector_1.push_back(-2.0);    vector_1.push_back(-1.5);
    vector_1.push_back(0.2);     vector_1.push_back(0.8);
    vector_1.push_back(1.3);
    EXPECT_TRUE( profugus::is_sorted(vector_1.begin(), vector_1.end()) );

    // Test a vector of non-sorted doubles
    Vec_Dbl vector_2;    // {0.2, 0.6, -1.2, 0.8, 1.0}
    vector_2.push_back(0.2);     vector_2.push_back(0.6);
    vector_2.push_back(-1.2);    vector_2.push_back(0.8);
    vector_2.push_back(1.0);
    EXPECT_TRUE( !profugus::is_sorted(vector_2.begin(), vector_2.end()) );

    // Test an unsorted deque of doubles
    std::deque<double> deque_1; // {0.6, 1.5, 0.9, 3.2, -1.0}
    deque_1.push_back(0.6);  deque_1.push_back(1.5);  deque_1.push_back(0.9);
    deque_1.push_back(3.2);  deque_1.push_back(-1.0);
    EXPECT_TRUE( !profugus::is_sorted(deque_1.begin(), deque_1.end()) );

    // Test a vector of sorted integers
    Vec_Int vector_3;   // {-3, -2, -1, 1, 5, 10}
    vector_3.push_back(-3);     vector_3.push_back(-2);
    vector_3.push_back(-1);     vector_3.push_back(1);
    vector_3.push_back(5);      vector_3.push_back(10);
    EXPECT_TRUE( profugus::is_sorted(vector_3.begin(), vector_3.end()) );

    // Test a vector of sorted strings
    Vec_String str_vector_1;
    str_vector_1.push_back("aaa");
    str_vector_1.push_back("bbb");
    str_vector_1.push_back("ccc");
    EXPECT_TRUE( profugus::is_sorted(str_vector_1.begin(),
                                     str_vector_1.end()) );

    // Test a vector of non-sorted strings
    Vec_String str_vector_2;
    str_vector_2.push_back("aaa");
    str_vector_2.push_back("bbb");
    str_vector_2.push_back("aaa");
    str_vector_2.push_back("ccc");
    EXPECT_TRUE( !profugus::is_sorted(str_vector_2.begin(),
                                      str_vector_2.end()) );

    // Test an edge case -- equal values
    Vec_Dbl vector_4(4, 4.0);
    EXPECT_FALSE( profugus::is_sorted(vector_4.begin(), vector_4.end(),
                                      std::less<double>()) );
}

//---------------------------------------------------------------------------//
TEST(IsSubset, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;
    typedef def::Vec_String Vec_String;

    // Make vec1 {0.0, 1.1, 2.2, 3.3, 4.4, 5.5}
    Vec_Dbl dbl_vec1(6, 0.0);
    for (int i = 0; i < dbl_vec1.size(); ++i)
    {
        dbl_vec1[i] = i * 1.1;
    }

    // Make vec2 {1.3, 3.3, 6.5}
    Vec_Dbl dbl_vec2(3);
    dbl_vec2[0] = 1.3; dbl_vec2[1] = 3.3; dbl_vec2[2] = 6.5;

    // Make vec3 {1.1, 3.3, 5.5}
    Vec_Dbl dbl_vec3(3);
    dbl_vec3[0] = 1.1; dbl_vec3[1] = 3.3; dbl_vec3[2] = 5.5;

    // vec2 is not a subset of vec1
    EXPECT_TRUE( !profugus::is_subset(dbl_vec1, dbl_vec2) );
    // vec3 is a subset of vec1
    EXPECT_TRUE( profugus::is_subset(dbl_vec1, dbl_vec3, 1.0e-6) );

    // Make Vec1 (superset of vec2)
    Vec_String str_vec1;
    str_vec1.push_back("aaa");
    str_vec1.push_back("bbb");
    str_vec1.push_back("ccc");
    str_vec1.push_back("ddd");
    str_vec1.push_back("eee");

    // Make Vec2 (subset of vec1)
    Vec_String str_vec2;
    str_vec2.push_back("aaa");
    str_vec2.push_back("ccc");
    str_vec2.push_back("eee");
    EXPECT_TRUE( profugus::is_subset(str_vec1, str_vec2) );

    // Make Vec3 (not a subset of vec1)
    Vec_String str_vec3;
    str_vec3.push_back("aaa");
    str_vec3.push_back("ccc");
    str_vec3.push_back("fff");
    EXPECT_TRUE( !profugus::is_subset(str_vec1, str_vec3) );
}

//---------------------------------------------------------------------------//
// Test is_unique
TEST(IsUnique, all)
{
    typedef def::Vec_Dbl    Vec_Dbl;
    typedef def::Vec_String Vec_Str;

    // Create a vector with only unique entries to test
    double dbl_test1[] = {1.0, 3.0, 5.0, 9.0, 11.0, 13.0, 16.0, 18.0, 21.0};
    Vec_Dbl dbl_test1_v(&dbl_test1[0], &dbl_test1[0]+9);
    EXPECT_TRUE( profugus::is_unique(dbl_test1_v.begin(), dbl_test1_v.end()) );

    // Create a vector with a repeat entry to test
    double dbl_test2[] = {1.0, 3.0, 5.0, 9.0, 11.0, 11.000001, 16, 18, 21};
    Vec_Dbl dbl_test2_v(&dbl_test2[0], &dbl_test2[0]+9);
    EXPECT_FALSE( profugus::is_unique(dbl_test2_v.begin(), dbl_test2_v.end(),
                                    0.001) );

    // Create a vector of strings with only unique entries to test
    Vec_Str str_test1_v;
    str_test1_v.push_back("aaa");
    str_test1_v.push_back("bbb");
    str_test1_v.push_back("ccc");
    str_test1_v.push_back("ddd");
    str_test1_v.push_back("eee");
    str_test1_v.push_back("fff");
    EXPECT_TRUE( profugus::is_unique(str_test1_v.begin(), str_test1_v.end()) );

    // Create a vector with a repeat entry to test
    Vec_Str str_test2_v;
    str_test2_v.push_back("aaa");
    str_test2_v.push_back("bbb");
    str_test2_v.push_back("ccc");
    str_test2_v.push_back("ddd");
    str_test2_v.push_back("aaa");
    str_test2_v.push_back("eee");
    str_test2_v.push_back("fff");
    EXPECT_FALSE( profugus::is_unique(str_test2_v.begin(), str_test2_v.end()) );
}

//---------------------------------------------------------------------------//
TEST(IsPositive, all)
{
    using profugus::is_positive;
    using profugus::constants::tiny;

    std::vector<double> v;
    EXPECT_TRUE(is_positive(v.begin(), v.end()));

    v.clear(); v.push_back(2.0);
    EXPECT_TRUE(is_positive(v.begin(), v.end()));

    v.clear(); v.push_back(-1.0);
    EXPECT_FALSE(is_positive(v.begin(), v.end()));

    v.clear(); v.push_back(2.0); v.push_back(-1.0);
    EXPECT_FALSE(is_positive(v.begin(), v.end()));

    v.clear(); v.push_back(-1.0); v.push_back(2.0);
    EXPECT_FALSE(is_positive(v.begin(), v.end()));

    v.clear(); v.push_back(0.0); v.push_back(2.0);
    EXPECT_FALSE(is_positive(v.begin(), v.end())) << "failure on 0.0";

    v.clear(); v.push_back(2.0); v.push_back(-0.0);
    EXPECT_FALSE(is_positive(v.begin(), v.end())) << "failure on -0.0";

    v.clear(); v.push_back(tiny);
    EXPECT_TRUE(is_positive(v.begin(), v.end())) << "failure on +precision";

    v.clear(); v.push_back(-1.0 * tiny);
    EXPECT_FALSE(is_positive(v.begin(), v.end())) << "failure on -precision";
}

//---------------------------------------------------------------------------//
TEST(IsNegative, all)
{
    using profugus::is_negative;
    using profugus::constants::tiny;

    std::vector<double> v;
    EXPECT_TRUE(is_negative(v.begin(), v.end()));

    v.clear(); v.push_back(2.0);
    EXPECT_FALSE(is_negative(v.begin(), v.end()));

    v.clear(); v.push_back(-1.0);
    EXPECT_TRUE(is_negative(v.begin(), v.end()));

    v.clear(); v.push_back(2.0); v.push_back(-1.0);
    EXPECT_FALSE(is_negative(v.begin(), v.end()));

    v.clear(); v.push_back(-1.0); v.push_back(2.0);
    EXPECT_FALSE(is_negative(v.begin(), v.end()));

    v.clear(); v.push_back(0.0); v.push_back(2.0);
    EXPECT_FALSE(is_negative(v.begin(), v.end())) << "failure on 0.0";

    v.clear(); v.push_back(2.0); v.push_back(-0.0);
    EXPECT_FALSE(is_negative(v.begin(), v.end())) << "failure on -0.0";

    v.clear(); v.push_back(tiny);
    EXPECT_FALSE(is_negative(v.begin(), v.end())) << "failure on +precision";

    v.clear(); v.push_back(-1.0 * tiny);
    EXPECT_TRUE(is_negative(v.begin(), v.end())) << "failure on -precision";
}

//---------------------------------------------------------------------------//
TEST(IsNonNegative, all)
{
    using profugus::is_non_negative;
    using profugus::constants::tiny;

    std::vector<double> v;
    EXPECT_TRUE(is_non_negative(v.begin(), v.end()));

    v.clear(); v.push_back(2.0);
    EXPECT_TRUE(is_non_negative(v.begin(), v.end()));

    v.clear(); v.push_back(-1.0);
    EXPECT_FALSE(is_non_negative(v.begin(), v.end()));

    v.clear(); v.push_back(2.0); v.push_back(-1.0);
    EXPECT_FALSE(is_non_negative(v.begin(), v.end()));

    v.clear(); v.push_back(-1.0); v.push_back(2.0);
    EXPECT_FALSE(is_non_negative(v.begin(), v.end()));

    v.clear(); v.push_back(0.0); v.push_back(2.0);
    EXPECT_TRUE(is_non_negative(v.begin(), v.end())) << "failure on 0.0";

    v.clear(); v.push_back(2.0); v.push_back(-0.0);
    EXPECT_TRUE(is_non_negative(v.begin(), v.end())) << "failure on -0.0";

    v.clear(); v.push_back(tiny);
    EXPECT_TRUE(is_non_negative(v.begin(), v.end())) << "failure on +precision";

    v.clear(); v.push_back(-1 * tiny);
    EXPECT_FALSE(is_non_negative(v.begin(), v.end()))
        << "failure on -precision";
}

//---------------------------------------------------------------------------//
//                 end of tstContainer_Props.cc
//---------------------------------------------------------------------------//
