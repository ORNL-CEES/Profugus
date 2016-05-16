//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/test/tstView_Field.cc
 * \author Thomas M. Evans et al
 * \date   Thu Oct 20 13:26:04 2011
 * \brief  View_Field unit test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../View_Field.hh"
#include "../View_Field_Struct.hh"

#include "Utils/gtest/utils_gtest.hh"

#include <vector>
#include <cmath>
#include <sstream>

#include "Utils/utils/Definitions.hh"
#include "Utils/comm/Timer.hh"

using profugus::make_view;
typedef profugus::View_Field<double>       View_Field_Dbl;
typedef profugus::const_View_Field<double> const_View_Field_Dbl;

using namespace std;
using def::Vec_Dbl;
using def::Vec_Int;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(View, basic)
{
    // Test a default View_Field_Dbl
    View_Field_Dbl empty_vf;
    EXPECT_EQ(empty_vf.begin().ptr(), nullptr);
    EXPECT_EQ(empty_vf.end().ptr(), nullptr);
    EXPECT_EQ(empty_vf.data(), nullptr);
    EXPECT_EQ(empty_vf.pbegin(), nullptr);
    EXPECT_EQ(empty_vf.pend(), nullptr);

    // >>> Test a View_Field_Dbl with pointers

    // Make a field
    vector<double> y(4);
    fill(y.begin(), y.end(), 0.57);
    y[2] = 0.62;

    // get a view
    {
        View_Field_Dbl v(&y[0], &y[3]+1);
        EXPECT_EQ(4, v.size());
        EXPECT_TRUE( !v.empty() );

        EXPECT_EQ(0.57, v[0]);
        EXPECT_EQ(0.57, v[1]);
        EXPECT_EQ(0.62, v[2]);
        EXPECT_EQ(0.57, v[3]);

        // change a value
        v[1] = 0.23;
    }

    // check field
    EXPECT_EQ(0.57, y[0]);
    EXPECT_EQ(0.23, y[1]);
    EXPECT_EQ(0.62, y[2]);
    EXPECT_EQ(0.57, y[3]);
    EXPECT_EQ(0.57, y.front());
    EXPECT_EQ(0.57, y.back());

    // Get another view
    {
        View_Field_Dbl v(&y[0], &y[3]+1);
        EXPECT_EQ(4, v.size());

        // Now, get the second half of y
        View_Field_Dbl x = v.slice(2, 1);

        EXPECT_EQ(2, x.size());
        EXPECT_EQ(0.62, x[0]);
        EXPECT_EQ(0.57, x[1]);
        EXPECT_EQ(0.62, x.front());
        EXPECT_EQ(0.57, x.back());

        // Now, get the first half of y
        View_Field_Dbl x2 = v.slice(1, 0, 2);

        EXPECT_EQ(2, x2.size());
        EXPECT_EQ(0.57, x2[0]);
        EXPECT_EQ(0.23, x2[1]);

        // Now, get index 1 through 3
        View_Field_Dbl z = v.general_slice(1, 4);

        EXPECT_EQ(3, z.size());
        EXPECT_EQ(0.23, z[0]);
        EXPECT_EQ(0.62, z[1]);
        EXPECT_EQ(0.57, z[2]);
    }
}

TEST(View, reverse_iterators)
{
    typedef profugus::View_Field<int>       View_Field_Int;
    typedef profugus::const_View_Field<int> const_View_Field_Int;

    const int orig_data[] = {1, 2, 3, 4, 5};
    Vec_Int mutable_data(orig_data, orig_data + 5);

    View_Field_Int view = make_view(mutable_data);
    EXPECT_VEC_EQ(orig_data, view);

    std::copy(orig_data, orig_data + 5, view.rbegin());

    const int reversed_data[] = {5, 4, 3, 2, 1};
    EXPECT_VEC_EQ(reversed_data, view);

    const_View_Field_Int cvf(reversed_data, reversed_data + 5);
    const Vec_Int temp(cvf.rbegin(), cvf.rend());
    EXPECT_VEC_EQ(orig_data, temp);
}

//---------------------------------------------------------------------------//

TEST(Const, View)
{
    // >>> Test a const_View_Field_Dbl with pointers

    // Test a default View_Field_Dbl
    const_View_Field_Dbl empty_vf;
    EXPECT_EQ(empty_vf.begin().ptr(), nullptr);
    EXPECT_EQ(empty_vf.end().ptr(), nullptr);
    EXPECT_EQ(empty_vf.data(), nullptr);

    // Make a field
    vector<double> y(6);
    fill(y.begin(), y.end(), 0.82);
    y[2] = 0.88;
    y[5] = 0.86;

    // get a view
    {
        const_View_Field_Dbl v(&y[0], &y[5]+1);
        EXPECT_EQ(6, v.size());
        EXPECT_TRUE( !v.empty() );

        EXPECT_EQ(0.82, v[0]);
        EXPECT_EQ(0.82, v[1]);
        EXPECT_EQ(0.88, v[2]);
        EXPECT_EQ(0.82, v[3]);
        EXPECT_EQ(0.82, v[4]);
        EXPECT_EQ(0.86, v[5]);
        EXPECT_EQ(0.82, v.front());
        EXPECT_EQ(0.86, v.back());

        // Now, get the middle values of y
        const_View_Field_Dbl x = v.slice(2, 1);

        EXPECT_EQ(2, x.size());
        EXPECT_EQ(0.88, x[0]);
        EXPECT_EQ(0.82, x[1]);
        EXPECT_EQ(0.88, x.front());
        EXPECT_EQ(0.82, x.back());

        // Now, get the last 2 values of y
        const_View_Field_Dbl x2 = v.slice(1, 4, 6);

        EXPECT_EQ(2, x2.size());
        EXPECT_EQ(0.82, x2[0]);
        EXPECT_EQ(0.86, x2[1]);

        // Now, get index 1 through 3
        const_View_Field_Dbl z = v.general_slice(1, 4);

        EXPECT_EQ(3, z.size());
        EXPECT_EQ(0.82, z[0]);
        EXPECT_EQ(0.88, z[1]);
        EXPECT_EQ(0.82, z[2]);
    }

    // check field
    EXPECT_EQ(0.82, y[0]);
    EXPECT_EQ(0.82, y[1]);
    EXPECT_EQ(0.88, y[2]);
    EXPECT_EQ(0.82, y[3]);
    EXPECT_EQ(0.82, y[4]);
    EXPECT_EQ(0.86, y[5]);
}

//---------------------------------------------------------------------------//

TEST(View, iterator_construction)
{
    typedef std::vector<int>               Vec_Int;
    typedef profugus::View_Field<int>       View_Field_Int;
    typedef profugus::const_View_Field<int> const_View_Field_Int;

    Vec_Int data{1,1,2,3,5,8};

    // Construct from pointers
    View_Field_Int view(data.data(), data.data() + data.size());

    // Construct subset
    View_Field_Int subview(view.begin() + 2, view.begin() + 4);
    EXPECT_VEC_EQ((Vec_Int{2,3}), subview);

    // Construct const view
    const_View_Field_Int cview(view.begin() + 3, view.end());
    EXPECT_VEC_EQ((Vec_Int{3,5,8}), cview);
}

//---------------------------------------------------------------------------//
TEST(Array, all)
{
    // Build an array
    int vals[] = {1, 3, 5, 7, 9, 11};

    // Make a view field of i's
    profugus::View_Field<int> view = make_view(vals);
    EXPECT_EQ(view.size(), 6);
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_EQ(vals[i], view[i]);
    }

    // Make the const view field of i's
    profugus::const_View_Field<int> const_view = make_view(vals);
    EXPECT_EQ(const_view.size(), 6);
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_EQ(vals[i], const_view[i]);
    }
}

//---------------------------------------------------------------------------//
TEST(Stride, all)
{
    // Build an array
    int vals[] = { 1,  3,  5,  7,  9,
                  11, 13, 15, 17, 19,
                  21, 23, 25, 27, 29,
                  31, 33, 35, 37, 39,
                  41, 43, 45, 47, 49};

    // Stride 1 view
    profugus::View_Field<int> view = make_view(vals);
    int counter = 1;
    for (int d : view)
    {
        EXPECT_EQ(counter, d);
        counter += 2;
    }

    // Stride 5 view
    profugus::View_Field<int> view2 = view.strided_slice(0, view.size(), 5);
    counter = 1;
    for (int d : view2)
    {
        EXPECT_EQ(counter, d);
        counter += 10;
    }

    // Stride 5 const view
    profugus::const_View_Field<int> const_view2(view2);
    counter = 1;
    for (int d : const_view2)
    {
        EXPECT_EQ(counter, d);
        counter += 10;
    }

    // Stride 5 view alternate
    profugus::View_Field<int> view3 = view.strided_slice(5, 25, 5);
    EXPECT_EQ(5, view3.stride());
    EXPECT_EQ(4, view3.size());
    EXPECT_VEC_EQ(std::vector<int>({11,21,31,41}),
                  std::vector<int>(view3.begin(), view3.end()));

    // Take a strided slice of a strided view
    profugus::View_Field<int> view4 = view3.strided_slice(1, 5, 2);
    EXPECT_EQ(5 * 2, view4.stride());
    EXPECT_VEC_EQ(std::vector<int>({21,41}),
                  std::vector<int>(view4.begin(), view4.end()));

    // Make a copy
    profugus::View_Field<int> view4_copy(view4);
    EXPECT_EQ(5 * 2, view4_copy.stride());
    EXPECT_VEC_EQ(view4, view4_copy);

    // Take a strided slice of a strided view (const)
    profugus::const_View_Field<int> view3c(view3);
    profugus::const_View_Field<int> view4c = view3c.strided_slice(1, 5, 2);
    EXPECT_VEC_EQ(view4, view4c);

    // Make a copy
    profugus::const_View_Field<int> view4c_copy(view4c);
    EXPECT_EQ(5 * 2, view4c_copy.stride());
    EXPECT_VEC_EQ(view4, view4c_copy);
}

//---------------------------------------------------------------------------//
TEST(Stride, more_copy_constructors)
{
    typedef profugus::const_View_Field<int> const_view;
    typedef profugus::View_Field<int>       view;

    // Build an array
    int vals[] = {1, 2, 3, 4};
    view vals_view = make_view(vals);

    // Mutable view
    view mv = vals_view.strided_slice(0, vals_view.size(), 2);
    const view& mv_cref(mv);

    // Const view
    const_view cv = vals_view.strided_slice(0, vals_view.size(), 2);
    const const_view& cv_cref(cv);

    ASSERT_EQ(2, mv.stride());
    ASSERT_EQ(2, mv_cref.stride());
    ASSERT_EQ(2, cv.stride());
    ASSERT_EQ(2, cv_cref.stride());

    // view -> view
    EXPECT_EQ(2, view(mv).stride());
    EXPECT_EQ(2, view(mv_cref).stride());

    // view -> const_view
    EXPECT_EQ(2, const_view(mv).stride());
    EXPECT_EQ(2, const_view(mv_cref).stride());

    // const_view -> const_view
    EXPECT_EQ(2, const_view(cv).stride());
    EXPECT_EQ(2, const_view(cv_cref).stride());
}

//---------------------------------------------------------------------------//
TEST(Copy, all)
{
    // Build an array
    int vals_1[] = { 1,  3,  5,  7,  9,
                     11, 13, 15, 17, 19,
                     21, 23, 25, 27, 29,
                     31, 33, 35, 37, 39,
                     41, 43, 45, 47, 49};

    // Build an array
    int vals_2[] = { 100, 103, 105, 107, 109,
                     110, 130, 115, 117, 119,
                     121, 123, 125, 127, 129,
                     131, 133, 135, 137, 139,
                     141, 143, 145, 147, 149};

    // Make vectors
    std::vector<int> vec_1(std::begin(vals_1), std::end(vals_1));
    std::vector<int> vec_2(std::begin(vals_2), std::end(vals_2));

    // View Field assign test (stride 1)
    {
        // Make view of vec_1 with stride 1
        profugus::View_Field<int> view_1 = make_view(vec_1);
        // Make a view of vec_2 with stride 1
        profugus::View_Field<int> view_2 = make_view(vec_2);
        // Copy from view_1 to view_2
        view_2.assign(view_1.begin(), view_1.end());
        // Check
        for (int i = 0; i < view_1.size(); ++i)
        {
            EXPECT_EQ(view_1[i], view_2[i]);
        }
    }
    // Vector assign test
    {
        // Make a view of vec_1 with stride 1
        profugus::View_Field<int> view_1 = make_view(vec_1);
        // Assign from vec_2
        view_1.assign(vec_2.begin(), vec_2.end());
        // Check
        for (int i = 0; i < view_1.size(); ++i)
        {
            EXPECT_EQ(vec_2[i], view_1[i]);
        }
    }
}

//---------------------------------------------------------------------------//
struct EnergyPoint
{
    float  a;
    float  b;
    double c;
};

TEST(StructView, manual)
{
    std::vector<EnergyPoint> energies = {
        {0.0f, 0.0f, 1.0},
        {0.2f, 0.2f, 2.0},
        {0.8f, 1.0f, 3.0},
    };

    // Create a view of the structs
    profugus::const_View_Field<EnergyPoint> eview = make_view(energies);
    EXPECT_FLOAT_EQ(0.0f, eview[0].a);
    EXPECT_FLOAT_EQ(0.2f, eview[1].a);
    EXPECT_FLOAT_EQ(0.8f, eview[2].a);

    // Try a view of the first element
    EXPECT_EQ(8, sizeof(double));
    EXPECT_EQ(16, sizeof(EnergyPoint));

    profugus::const_View_Field<double> cview(
            &(eview.data()->c),
            &((eview.data() + eview.size())->c),
            sizeof(EnergyPoint) / sizeof(double));

    double expected_cdata[] = {1.0, 2.0, 3.0};

    EXPECT_VEC_EQ(expected_cdata, cview);
}

TEST(StructView, energy_stride_member)
{
    std::vector<EnergyPoint> energies = {
        {0.0f, 0.0f, 1.0},
        {0.2f, 0.2f, 2.0},
        {0.8f, 1.0f, 3.0},
    };

    // Try a view of the first element
    EXPECT_EQ(4, sizeof(float));
    EXPECT_EQ(16, sizeof(EnergyPoint));

    auto aview = PROFUGUS_MAKE_STRUCT_VIEW(energies, a);
    auto bview = PROFUGUS_MAKE_STRUCT_VIEW(energies, b);
    auto cview = PROFUGUS_MAKE_STRUCT_VIEW(make_view(energies), c);

    EXPECT_VEC_EQ((std::vector<float>{.0f, .2f, .8f}), aview);
    EXPECT_VEC_EQ((std::vector<float>{.0f, .2f, 1.0f}), bview);
    EXPECT_VEC_EQ((Vec_Dbl{1, 2, 3}), cview);
}

struct CosinePoint
{
    float a;
    float b;
    float c;
};

TEST(StructView, cosine_stride)
{
    std::vector<CosinePoint> cosines = {
        {0.0f, 0.0f, 1.0f},
        {0.2f, 0.2f, 2.0f},
        {0.8f, 1.0f, 3.0f},
    };

    // Try a view of the second element
    EXPECT_EQ(4, sizeof(float));
    EXPECT_EQ(12, sizeof(CosinePoint));

    auto bview = PROFUGUS_MAKE_STRUCT_VIEW(cosines, b);

    float expected_values[] = {0.0, 0.2, 1.0};

    EXPECT_VEC_EQ(expected_values, bview);
}

struct SmallStruct
{
    char a;
    char b;
    char c;
};

// This tests possible alignment issues with small structures.
TEST(StructView, char_stride)
{
    std::vector<SmallStruct> chars = {
        {'a', 'b', 'c'},
        {'u', 'v', 'w'},
        {'x', 'y', 'z'},
    };

    // Try a view of the third element
    EXPECT_EQ(1, sizeof(char));
    EXPECT_EQ(3, sizeof(SmallStruct));
    EXPECT_EQ(1, alignof(SmallStruct));

    auto cview = PROFUGUS_MAKE_STRUCT_VIEW(chars, c);

    char expected_values[] = {'c', 'w', 'z'};
    EXPECT_VEC_EQ(expected_values, cview);
}

struct WackyStruct
{
    char   a;
    int    b;
    double c;
};

// This tests other possible alignment issues?
TEST(StructView, wacky)
{
    std::vector<WackyStruct> vals = {
        {'a', 3,   100.},
        {'u', 5,  1000.},
        {'x', 8, 10000.},
    };

    // Padding gets built into struct size
    EXPECT_EQ(1 + 3 + sizeof(int) + sizeof(double), sizeof(WackyStruct));
    EXPECT_EQ(alignof(double), alignof(WackyStruct));

    // Try a view of the first element
    auto aview = PROFUGUS_MAKE_STRUCT_VIEW(vals, a);

    char expected_avalues[] = {'a', 'u', 'x'};
    EXPECT_VEC_EQ(expected_avalues, aview);

    // Try a view of the second element
    auto bview = PROFUGUS_MAKE_STRUCT_VIEW(vals, b);

    int expected_bvalues[] = {3, 5, 8};
    EXPECT_VEC_EQ(expected_bvalues, bview);

    // Try a view of the third element
    auto cview = PROFUGUS_MAKE_STRUCT_VIEW(vals, c);

    double expected_cvalues[] = {100, 1000, 10000};
    EXPECT_VEC_EQ(expected_cvalues, cview);
}

//---------------------------------------------------------------------------//

TEST(RangeView, range)
{
    Vec_Int vars{1, 2, 3, 6};
    auto view = profugus::make_view(vars);

    for (auto& v : view.fast_range())
    {
        v += 3;
    }

    EXPECT_VEC_EQ((Vec_Int{4, 5, 6, 9}), vars);

    // Check iterating on an empty field
    profugus::View_Field<int> empty_view;
    int ctr = 0;
    for (auto v : empty_view)
    {
        ++ctr;
    }
    EXPECT_EQ(0, ctr);

#ifdef REQUIRE_ON
    profugus::View_Field<int> strided_view(vars.data(),
                                          vars.data() + vars.size(),
                                          2);
    EXPECT_THROW(strided_view.fast_range(), profugus::assertion);
#endif
}

TEST(RangeView, piter)
{
    Vec_Int vars{1, 2, 3, 6};
    auto view = profugus::make_view(vars);

    EXPECT_EQ(view.data(), view.pbegin());
    EXPECT_EQ(view.data(), view.cpbegin());
    EXPECT_EQ(view.data() + view.size(), view.pend());
    EXPECT_EQ(view.data() + view.size(), view.cpend());

#ifdef REQUIRE_ON
    profugus::View_Field<int> strided_view(vars.data(),
                                          vars.data() + vars.size(),
                                          2);
    EXPECT_THROW(strided_view.pbegin(), profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
TEST(Timing, DISABLED_stride_1)
{
    // Create a vector of source data
    Vec_Dbl source(500000, 0.0);
    for (int i = 0; i < 500000; ++i)
    {
        source[i] = static_cast<double>(i);
    }

    // Create a view of source data
    profugus::const_View_Field<double> source_view = make_view(source);

    // Create another vector of data
    Vec_Dbl target(500000, 0.0);

    // Create a view of the target
    profugus::View_Field<double> target_view = make_view(target);

    Vec_Dbl times;
    profugus::Timer timer;

    for (int i = 0; i < 10; ++i)
    {

        timer.start();

        // Copy 1000 times
        for (int j = 0; j < 1000; ++j)
        {
            // Copy the data
            target_view.assign(source_view.begin(), source_view.end());
            /*
            std::copy(source_view.begin(), source_view.end(),
                      target_view.begin());
            */

            // Clear the target
            std::fill(target_view.begin(), target_view.end(), 0.0);
        }

        timer.stop();
        times.push_back(timer.wall_clock());
    }

    for (double d : times)
    {
        std::cout << "Time: " << d << std::endl;
    }
}

//---------------------------------------------------------------------------//
//                        end of tstView_Field.cc
//---------------------------------------------------------------------------//
