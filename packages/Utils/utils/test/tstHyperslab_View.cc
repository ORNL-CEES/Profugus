//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/test/tstHyperslab_View.cc
 * \author Seth R Johnson
 * \date   Mon Dec 08 16:14:17 2014
 * \brief  Test for Hyperslab_View
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Hyperslab_View.hh"

#include "Utils/gtest/utils_gtest.hh"

using profugus::Hyperslab_View;
using profugus::View_Field;
using profugus::const_Hyperslab_View;
using profugus::const_View_Field;

//---------------------------------------------------------------------------//
// HYPERSLAB TESTS
//---------------------------------------------------------------------------//

TEST(ConstHyperslabViewTest, two_dimension)
{
    typedef std::vector<int>              vec_t;
    typedef const_Hyperslab_View<int, 2>  const_hyperslab_t;
    typedef const_hyperslab_t::index_type index_t;

    vec_t temp = {1,2,3,4, 11,12,13,14, 21,22,23,24};

    index_t dims(3, 4); // 3 (row) major, 4 (column) minor
    const_hyperslab_t slab(profugus::make_view(temp), dims);

    EXPECT_EQ(slab[index_t(2, 3)], 24);
    // Slightly more expensive but easier to write method
    EXPECT_EQ(slab[2][3], 24);

#ifdef REQUIRE_ON
    EXPECT_THROW(slab[3][0], profugus::assertion);
    EXPECT_THROW(slab[0][4], profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
TEST(HyperslabViewTest, two_dimension)
{
    typedef std::vector<int>             vec_t;
    typedef const_Hyperslab_View<int, 2> const_hyperslab_t;
    typedef Hyperslab_View<int, 2>       hyperslab_t;
    typedef hyperslab_t::index_type      index_t;

    vec_t temp = {1,2,3,4, 11,12,13,14, 21,22,23,24};

    index_t dims(3, 4); // 3 (row) major, 4 (column) minor
    hyperslab_t slab(profugus::make_view(temp), dims);

    // Modify the data
    slab[index_t(2,3)] = -100;
    EXPECT_EQ(-100, temp[11]);
    slab[0][1] = -1;
    EXPECT_EQ(-1, temp[1]);

    // Modify as a view
    std::fill(slab[0].begin(), slab[0].end(), 999);
    EXPECT_EQ(999, temp[3]);
    EXPECT_EQ(11,  temp[4]);

    // Check implicit conversion to const slab
    auto get_val = [](const const_hyperslab_t &s){ return s[1][3]; };
    EXPECT_EQ(14, get_val(slab));
}

//---------------------------------------------------------------------------//
TEST(HyperslabViewTest, assignment)
{
    typedef std::vector<int>        vec_t;
    typedef Hyperslab_View<int, 2>  hyperslab_t;
    typedef hyperslab_t::index_type index_t;

    vec_t temp = {1,2,3,4, 11,12,13,14, 21,22,23,24};

    hyperslab_t slab;
    EXPECT_EQ(0, slab.size());

    index_t dims(3, 4); // 3 (row) major, 4 (column) minor
    hyperslab_t slab2(profugus::make_view(temp), dims);

    // Assignment
    slab = slab2;
    EXPECT_EQ(3*4, slab.size());
    EXPECT_EQ(dims, slab.dims());
    slab[index_t(1,1)] = 100;
    EXPECT_EQ(100, slab2[index_t(1,1)]);
}

//---------------------------------------------------------------------------//
TEST(HyperslabViewTest, strides)
{
    typedef std::vector<int>             vec_t;
    typedef const_Hyperslab_View<int, 2> const_hyperslab_t;
    typedef Hyperslab_View<int, 2>       hyperslab_t;
    typedef hyperslab_t::index_type      index_t;
    typedef View_Field<int>              view;

    vec_t data(4 * 5);
    {
        auto it = data.begin();
        for (int i = 0; i != 4; ++i)
        {
            for (int j = 0; j != 5; ++j)
            {
                *it++ = 10*(i + 1) + (j + 1);
            }
        }
        CHECK(it == data.end());
    }

    view data_view(data.data() + 1,
                   data.data() + data.size() + 1,
                   5);
    const int expected[] = {12, 22, 32, 42};
    EXPECT_EQ(5, data_view.stride());
    EXPECT_VEC_EQ(expected, data_view);

    // View hyperslab
    hyperslab_t mh(data_view, index_t(4, 1));
    const hyperslab_t &mh_cref(mh);

    // Const view hyperslab
    const_hyperslab_t ch(data_view, index_t(4, 1));
    const const_hyperslab_t &ch_cref(ch);

    ASSERT_EQ(5, mh.stride());
    ASSERT_EQ(5, mh_cref.stride());
    ASSERT_EQ(5, ch.stride());
    ASSERT_EQ(5, ch_cref.stride());

    // view -> view
    EXPECT_EQ(5, hyperslab_t(mh).stride());
    EXPECT_EQ(5, hyperslab_t(mh_cref).stride());

    // view -> const_view
    EXPECT_EQ(5, const_hyperslab_t(mh).stride());
    EXPECT_EQ(5, const_hyperslab_t(mh_cref).stride());

    // const  view -> const_view
    EXPECT_EQ(5, const_hyperslab_t(ch).stride());
    EXPECT_EQ(5, const_hyperslab_t(ch_cref).stride());

    // Minor slice on view
    auto oned_slab = mh.minor_slice(0);
    EXPECT_EQ(5, oned_slab.stride());
    EXPECT_VEC_EQ(expected, oned_slab);

    // Minor slice on const view
    auto oned_cslab = ch.minor_slice(0);
    EXPECT_EQ(5, oned_cslab.stride());
    EXPECT_VEC_EQ(expected, oned_cslab);
}

//---------------------------------------------------------------------------//
TEST(HyperslabViewTest, minor_slice)
{
    typedef std::vector<int>             vec_t;
    typedef View_Field<int>              view;
    typedef const_View_Field<int>        const_view;
    typedef const_Hyperslab_View<int, 4> const_hyperslab_t;
    typedef Hyperslab_View<int, 4>       hyperslab_t;
    typedef hyperslab_t::index_type      index_t;

    // Create and assign data
    vec_t data(2*3*4*5);
    index_t dims(2, 3, 4, 5);
    {
        index_t i;
        auto it = data.begin();
        for (i[0] = 0; i[0] != dims[0]; ++i[0])
        {
            for (i[1] = 0; i[1] != dims[1]; ++i[1])
            {
                for (i[2] = 0; i[2] != dims[2]; ++i[2])
                {
                    for (i[3] = 0; i[3] != dims[3]; ++i[3])
                    {
                        *it++ = (1000 * (i[0]+1)
                                + 100 * (i[1]+1)
                                +  10 * (i[2]+1)
                                +   1 * (i[3]+1));
                    }
                }
            }
        }
        CHECK(it == data.end());
    }

    // Create slab
    hyperslab_t slab(profugus::make_view(data), dims);

    // Take a slice on the minor axis
    {
        typedef const_Hyperslab_View<int, 3> const_subslab_t;
        typedef const_subslab_t::index_type  subindex_t;
        const_subslab_t subslab = const_hyperslab_t(slab).minor_slice(3);

        subindex_t dims = subslab.dims();
        EXPECT_VEC_EQ(subindex_t(2, 3, 4), dims);

        subindex_t i;
        auto it = subslab.begin();
        EXPECT_EQ(5, it.stride());
        for (i[0] = 0; i[0] != dims[0]; ++i[0])
        {
            for (i[1] = 0; i[1] != dims[1]; ++i[1])
            {
                for (i[2] = 0; i[2] != dims[2]; ++i[2])
                {
                    EXPECT_EQ(*it++, (1000 * (i[0] + 1)
                                     + 100 * (i[1] + 1)
                                     +  10 * (i[2] + 1)
                                     +   1 * (3    + 1)));
                }
            }
        }
        EXPECT_TRUE(it == subslab.end());

        auto subsubslab = subslab.minor_slice(2);
        it = subsubslab.begin();
        EXPECT_EQ(5 * 4, it.stride());
        for (i[0] = 0; i[0] != dims[0]; ++i[0])
        {
            for (i[1] = 0; i[1] != dims[1]; ++i[1])
            {
                EXPECT_EQ(*it++, (1000 * (i[0] + 1)
                                 + 100 * (i[1] + 1)
                                 +  10 * (2    + 1)
                                 +   1 * (3    + 1)));
            }
        }
        EXPECT_TRUE(it == subsubslab.end());
    }

    // Multiple subslices with mutable slab
    const int expected[] = {1123, 1223, 1323, 2123, 2223, 2323};
    {
        auto temp = slab.minor_slice(3 - 1).minor_slice(2 - 1);
        EXPECT_EQ(5 * 4, temp.begin().stride());
        EXPECT_EQ(5 * 1 + 2,  temp.begin().ptr()
                            - slab.begin().ptr());
        EXPECT_VEC_EQ(expected, temp);

        // Check that view creation from mutable subslabs works
        view subview(temp);
        EXPECT_EQ(5 * 4, subview.begin().stride());
        EXPECT_EQ(5 * 1 + 2,  subview.begin().ptr()
                            - slab.begin().ptr());
    }
    {
        const hyperslab_t& const_slab(slab);
        auto temp = const_slab.minor_slice(3 - 1).minor_slice(2 - 1);
        EXPECT_EQ(5 * 4, temp.begin().stride());
        EXPECT_EQ(5 * 1 + 2,  temp.begin().ptr()
                            - slab.begin().ptr());
        EXPECT_VEC_EQ(expected, temp);
    }
    {
        const_hyperslab_t const_slab(slab);
        auto temp = const_slab.minor_slice(3 - 1).minor_slice(2 - 1);
        EXPECT_EQ(5 * 4, temp.begin().stride());
        EXPECT_EQ(5 * 1 + 2,  temp.begin().ptr()
                            - const_slab.begin().ptr());
        EXPECT_VEC_EQ(expected, temp);

        // Check that view creation from subslabs works
        const_view subview(temp);
        EXPECT_EQ(5 * 4, subview.begin().stride());
        EXPECT_EQ(5 * 1 + 2,  subview.begin().ptr()
                            - const_slab.begin().ptr());
    }
}

//---------------------------------------------------------------------------//
//                 end of tstHyperslab_View.cc
//---------------------------------------------------------------------------//
