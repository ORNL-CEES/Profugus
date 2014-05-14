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

#include <memory>
#include <vector>

//---------------------------------------------------------------------------//
// Test helpers
//---------------------------------------------------------------------------//

class Foo
{
  public:
    enum Foos
    {
        A, B, C
    };

  private:
    std::vector<int> d_data;

  public:
    Foo() : d_data(10, 15) {}
    const std::vector<int>& data() const { return d_data; }
    std::vector<int>& data() { return d_data; }

    Foos foo_type() const;

    std::shared_ptr<std::vector<int>> state() const
    {
        auto x = std::make_shared<decltype(d_data)>(d_data);
        return x;
    }
};

auto Foo::foo_type() const -> Foos
{
    return A;
}

template<class T>
auto change_data(const T &f) -> decltype(f.state())
{
    auto t = f.state();
    for (auto &i : *t)
    {
        i += 1;
    }
    return t;
}

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

    // decltype
    auto t = f.foo_type();
    EXPECT_EQ(Foo::Foos::A, t);
}

//---------------------------------------------------------------------------//

TEST(DecltypeTest, define)
{
    int y = 10;
    decltype(y) x = y;
    EXPECT_EQ(10, x);

    // this is the same
    auto z = y;
    EXPECT_EQ(10, z);
}

//---------------------------------------------------------------------------//

TEST(DecltypeTest, decltype)
{
    Foo f;

    auto d = change_data(f);

    for (auto i : *d)
    {
        EXPECT_EQ(16, i);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstAuto.cc
//---------------------------------------------------------------------------//
