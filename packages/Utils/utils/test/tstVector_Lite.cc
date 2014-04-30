//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstVector_Lite.cc
 * \author Thomas M. Evans
 * \date   Thu Jan  3 11:52:10 2008
 * \brief  Vector_Lite test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>

#include "../Vector_Lite.hh"

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

    {
        cout << "fill in from C array" << endl;
        double v1[5];
        for (size_t i=0; i<5; i++) {v1[i] = 1.*i;}
        Vector_Lite<double, 5> v2; v2.fill(v1);
        for (size_t i=0; i<5; i++) {
            EXPECT_EQ(v2[i], v1[i]);
        }
    }

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
        EXPECT_TRUE(profugus::soft_equiv(z.begin(), z.end(),
                                       ans.begin(), ans.end()));

        cout << "operator/=" << endl;
        z /= y;
        cout << " z = " << z << endl;
        EXPECT_TRUE(profugus::soft_equiv(z.begin(), z.end(),
                                       x.begin(), x.end()));
    }

    {
        cout << "operator+=" << endl;
        Vector_Lite<double, m> z(x);
        Vector_Lite<double, m> ans(c1+y0, c1+y1, c1+y2, c1+y3, c1+y4);
        z += y;
        cout << " z = " << z << endl;
        EXPECT_TRUE(profugus::soft_equiv(z.begin(), z.end(),
                                       ans.begin(), ans.end()));

        cout << "operator-=" << endl;
        z -= y;
        cout << " z = " << z << endl;
        EXPECT_TRUE(profugus::soft_equiv(z.begin(), z.end(),
                                       x.begin(), x.end()));
    }

    {
        cout << "unary-" << endl;
        Vector_Lite<double, m> z;
        Vector_Lite<double, m> ans(-c1);
        z = -x;
        cout << " z = " << z << endl;
        EXPECT_TRUE(profugus::soft_equiv(z.begin(), z.end(),
                                      ans.begin(), ans.end()));
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
        xNest(1) = y;
        cout << xNest << endl;
        EXPECT_EQ(x, xNest(0));
        EXPECT_EQ(y, xNest(1));
        EXPECT_EQ(x, xNest(2));
    }

    int i = -1;
    cout << "Negative bounds check x(" << i << ")\n";
    try {
        x(i);
    }
    catch ( profugus::assertion &a ) {
        EXPECT_TRUE(1);
    }
    catch (...) {
        cout << "Unknown error thrown.\n";
        EXPECT_TRUE(0);
    }

    Vector_Lite<double, m>::size_type iu(i);
    cout << "Negative bounds check test, unsigned x(" << iu << ")\n";
    try {
        x(iu);
    }
    catch ( profugus::assertion &a ) {
        EXPECT_TRUE(1);
    }
    catch (...) {
        cout << "Unknown error thrown.\n";
        EXPECT_TRUE(0);
    }

    i = x.size();
#ifdef REQUIRE_ON
    EXPECT_THROW(x(i), profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
//                        end of tstVector_Lite.cc
//---------------------------------------------------------------------------//
