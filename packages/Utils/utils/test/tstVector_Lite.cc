//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/test/tstVector_Lite.cc
 * \author Thomas M. Evans
 * \date   Thu Jan  3 11:52:10 2008
 * \brief  Vector_Lite test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Vector_Lite.hh"

#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <functional>

#include "Utils/gtest/utils_gtest.hh"
#include "Utils/utils/Hash_Functions.hh"

using namespace std;
using namespace profugus;

using profugus::Vector_Lite;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(VectorLite, basic)
{
    const int m = 5;

    cout << "constructor from scalar" << endl;
    Vector_Lite<double, m> x(0.0);
    EXPECT_EQ(m, std::count(x.begin(), x.end(), 0.0));

    cout << "assignment from another Vector_Lite" << endl;
    Vector_Lite<int, 3> ix(0, 1, 2);
    Vector_Lite<int, 3> iy(5, 6, 7);
    iy    = ix;
    ix[1] = 4;
    {
        EXPECT_EQ(0, ix[0]);
        EXPECT_EQ(4, ix[1]);
        EXPECT_EQ(2, ix[2]);
        EXPECT_EQ(0, iy[0]);
        EXPECT_EQ(1, iy[1]);
        EXPECT_EQ(2, iy[2]);
    }

    cout << "assignment to scalar" << endl;
    double c1 = 3.0;
    x = c1;
    cout << "x = " << x << endl;
    EXPECT_EQ(m, std::count(x.begin(), x.end(), c1));

    cout << "operator==" << endl;
    EXPECT_EQ(x, x);

    {
        cout << "copy constructor" << endl;
        Vector_Lite<double, m> xCopy(x);
        EXPECT_EQ(xCopy, x);
    }

    cout << "operator+=, scalar" << endl;
    double dc1 = 2.3;
    c1 += dc1;
    x += dc1;
    cout << " x = " << x << endl;
    EXPECT_EQ(m, std::count(x.begin(), x.end(), c1));

    cout << "operator-=, scalar" << endl;
    c1 -= dc1;
    x -= dc1;
    cout << " x = " << x << endl;
    EXPECT_EQ(m, std::count(x.begin(), x.end(), c1));

    cout << "operator*=, scalar" << endl;
    c1 *= dc1;
    x *= dc1;
    cout << " x = " << x << endl;
    EXPECT_EQ(m, std::count(x.begin(), x.end(), c1));

    cout << "operator/=, scalar" << endl;
    c1 /= dc1;
    x /= dc1;
    cout << " x = " << x << endl;
    EXPECT_EQ(m, std::count(x.begin(), x.end(), c1));

    double y0 = 2.0;
    double y1 = 1.0;
    double y2 = 0.3;
    double y3 = 0.2;
    double y4 = 62.7;
    Vector_Lite<double, m> y(y0, y1, y2, y3, y4);

    {
        cout << "operator*=" << endl;
        Vector_Lite<double, m> z(x);
        Vector_Lite<double, m> ans(c1*y0, c1*y1, c1*y2, c1*y3, c1*y4);
        z *= y;
        cout << " z = " << z << endl;
        EXPECT_VEC_SOFT_EQ(ans, z);

        cout << "operator/=" << endl;
        z /= y;
        cout << " z = " << z << endl;
        EXPECT_VEC_SOFT_EQ(x, z);
    }

    {
        cout << "operator+=" << endl;
        Vector_Lite<double, m> z(x);
        Vector_Lite<double, m> ans(c1+y0, c1+y1, c1+y2, c1+y3, c1+y4);
        z += y;
        EXPECT_VEC_SOFT_EQ(ans, z);

        cout << "operator-=" << endl;
        z -= y;
        EXPECT_VEC_SOFT_EQ(x, z);
    }

    {
        cout << "unary-" << endl;
        Vector_Lite<double, m> z;
        Vector_Lite<double, m> ans(-c1);
        z = -x;
        EXPECT_VEC_SOFT_EQ(ans, z);
    }

    {
        cout << "Inner product" << endl;
        Vector_Lite<double, 2> x1(1.0, 2.0);
        Vector_Lite<double, 2> x2(4.0, 6.0);
        EXPECT_EQ(16.0, profugus::inner_product(x1, x2));
    }

    {
        cout << "Nested Vector_Lites ";
        Vector_Lite<Vector_Lite<double, m>, 3> xNest(x);
        xNest[1] = y;
        cout << xNest << endl;
        EXPECT_EQ(x, xNest[0]);
        EXPECT_EQ(y, xNest[1]);
        EXPECT_EQ(x, xNest[2]);
    }

#ifdef REQUIRE_ON
    EXPECT_THROW(x[x.size()], profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//

TEST(VectorLite, comparison)
{
    typedef Vector_Lite<int, 3> vec_t;

    // lex comparison
    EXPECT_TRUE( vec_t(1, 2, 3) < vec_t(2, 1, 1));
    EXPECT_FALSE(vec_t(1, 2, 3) < vec_t(1, 2, 3));
    EXPECT_FALSE(vec_t(1, 2, 3) < vec_t(0, 2, 3));

    EXPECT_TRUE( vec_t(2, 1, 1) > vec_t(1, 2, 3));
    EXPECT_FALSE(vec_t(1, 2, 3) > vec_t(1, 2, 3));
    EXPECT_FALSE(vec_t(0, 2, 3) > vec_t(1, 2, 3));

    EXPECT_FALSE(vec_t(1, 2, 3) <= vec_t(1, 1, 1));
    EXPECT_TRUE( vec_t(1, 2, 3) <= vec_t(2, 1, 1));
    EXPECT_TRUE( vec_t(1, 2, 3) <= vec_t(1, 2, 3));

    EXPECT_FALSE(vec_t(1, 1, 1) >= vec_t(1, 2, 3));
    EXPECT_TRUE( vec_t(2, 1, 1) >= vec_t(1, 2, 3));
    EXPECT_TRUE( vec_t(1, 2, 3) >= vec_t(1, 2, 3));

    // "all" comparison
    EXPECT_TRUE( vec_t(1, 2, 3).all_lt(vec_t(2, 3, 4)));
    EXPECT_FALSE(vec_t(1, 2, 3).all_lt(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(1, 2, 3).all_lt(vec_t(0, 2, 3)));
    EXPECT_FALSE(vec_t(1, 2, 3).all_lt(vec_t(2, 1, 1)));

    EXPECT_TRUE( vec_t(2, 3, 4).all_gt(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(1, 2, 3).all_gt(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(0, 2, 3).all_gt(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(2, 1, 1).all_gt(vec_t(1, 2, 3)));

    EXPECT_TRUE( vec_t(1, 2, 3).all_le(vec_t(2, 3, 4)));
    EXPECT_TRUE( vec_t(1, 2, 3).all_le(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(1, 2, 3).all_le(vec_t(0, 2, 3)));
    EXPECT_FALSE(vec_t(1, 2, 3).all_le(vec_t(2, 1, 1)));

    EXPECT_TRUE( vec_t(2, 3, 4).all_ge(vec_t(1, 2, 3)));
    EXPECT_TRUE( vec_t(1, 2, 3).all_ge(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(0, 2, 3).all_ge(vec_t(1, 2, 3)));
    EXPECT_FALSE(vec_t(2, 1, 1).all_ge(vec_t(1, 2, 3)));
}

//---------------------------------------------------------------------------//

TEST(VectorLite, reverse_iterators)
{
    typedef std::vector<int> Vec_Int;

    Vector_Lite<int, 5> orig_data(1, 2, 3, 4, 5);
    Vector_Lite<int, 5> reversed_data(5, 4, 3, 2, 1);

    Vec_Int mutable_data(orig_data.rbegin(), orig_data.rend());
    EXPECT_VEC_EQ(reversed_data, mutable_data);
}

//---------------------------------------------------------------------------//

TEST(VectorLite, initializer_lists)
{
    Vector_Lite<int, 1> a = {1};
    Vector_Lite<int, 2> b = {1, 2};
    Vector_Lite<int, 3> c = {1, 2, 3};
    Vector_Lite<int, 4> d = {1, 2, 3, 4};
    Vector_Lite<int, 5> e = {1, 2, 3, 4, 5};
    Vector_Lite<int, 6> f = {1, 2, 3, 4, 5, 6};
    Vector_Lite<int, 7> g = {1, 2, 3};

    EXPECT_EQ(1, a[0]);

    EXPECT_EQ(1, b[0]);
    EXPECT_EQ(2, b[1]);

    EXPECT_EQ(1, c[0]);
    EXPECT_EQ(2, c[1]);
    EXPECT_EQ(3, c[2]);

    EXPECT_EQ(1, d[0]);
    EXPECT_EQ(2, d[1]);
    EXPECT_EQ(3, d[2]);
    EXPECT_EQ(4, d[3]);

    EXPECT_EQ(1, e[0]);
    EXPECT_EQ(2, e[1]);
    EXPECT_EQ(3, e[2]);
    EXPECT_EQ(4, e[3]);
    EXPECT_EQ(5, e[4]);

    EXPECT_EQ(1, f[0]);
    EXPECT_EQ(2, f[1]);
    EXPECT_EQ(3, f[2]);
    EXPECT_EQ(4, f[3]);
    EXPECT_EQ(5, f[4]);
    EXPECT_EQ(6, f[5]);

    EXPECT_EQ(1, g[0]);
    EXPECT_EQ(2, g[1]);
    EXPECT_EQ(3, g[2]);
    EXPECT_EQ(0, g[3]);
    EXPECT_EQ(0, g[4]);
    EXPECT_EQ(0, g[5]);
    EXPECT_EQ(0, g[6]);

#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        Vector_Lite<int, 2> h{1, 2, 3};
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}

//---------------------------------------------------------------------------//
//                        end of tstVector_Lite.cc
//---------------------------------------------------------------------------//
