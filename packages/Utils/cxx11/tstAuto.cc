//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cxx11/tstAuto.cc
 * \author Thomas M. Evans
 * \date   Fri Mar 14 17:02:20 2014
 * \brief  CXX-11 Auto test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Foo
{
  private:
    std::vector<int> d_data;

  public:
    Foo() : d_data(10, 15) {}
    const std::vector<int>& data() const { return d_data; }
    std::vector<int>& data() { return d_data; }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(AutoTest, iterators)
{
    int i = 1;
    std::vector<int> x(5);

    for (auto itr = x.begin(); itr != x.end(); ++itr)
    {
        *itr = i;
        ++i;
    }

    int ref[] = {1, 2, 3, 4, 5};
    EXPECT_VEC_EQ(ref, x);
}

//---------------------------------------------------------------------------//

TEST(AutoTest, return_type)
{
    Foo f;

    // constant reference
    const auto &d = f.data();

    int ref[] = {15, 15, 15, 15, 15, 15, 15, 15, 15, 15};

    EXPECT_VEC_EQ(ref, d);

    // mutable reference
    auto &m = f.data();
    m[3]    = 12;

    int refm[] = {15, 15, 15, 12, 15, 15, 15, 15, 15, 15};

    EXPECT_VEC_EQ(refm, f.data());

    // copy
    auto c = f.data();
    c[3]   = 11;

    int refc[] = {15, 15, 15, 11, 15, 15, 15, 15, 15, 15};

    EXPECT_VEC_EQ(refm, f.data());
    EXPECT_VEC_EQ(refc, c);
}

//---------------------------------------------------------------------------//
//                 end of tstAuto.cc
//---------------------------------------------------------------------------//
