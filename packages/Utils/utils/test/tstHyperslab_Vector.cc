//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/test/tstHyperslab_Vector.cc
 * \author Seth R Johnson
 * \date   Mon Aug 10 16:33:13 2015
 * \brief  Hyperslab_Vector class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Hyperslab_Vector.hh"

#include "Utils/gtest/utils_gtest.hh"

using profugus::Hyperslab_Vector;
using profugus::make_view;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class HyperslabVectorTest : public ::testing::Test
{
  protected:
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HyperslabVectorTest, subtest_description)
{
    typedef Hyperslab_Vector<int, 2> Hypervec_t;
    typedef Hypervec_t::const_View_t const_View_t;
    typedef std::vector<int>         Vec_t;
    typedef Hypervec_t::index_type   index_type;

    Hypervec_t h;
    EXPECT_EQ(0, h.size());
    EXPECT_TRUE(h.empty());

    int data[] = { 1,  2,  3,
                  11, 12, 13};

    h = Hypervec_t(const_View_t(make_view(data), index_type(2, 3)));

    EXPECT_EQ(6, h.size());
    EXPECT_EQ(index_type(2, 3), h.dims());

    EXPECT_VEC_EQ((Vec_t{11, 12, 13}), h.cview()[1]);
    EXPECT_VEC_EQ((Vec_t{3, 13}), h.cview().minor_slice(2));

    h.view()[index_type(1,1)] = 99;
    EXPECT_VEC_EQ((Vec_t{11, 99, 13}), h.view()[1]);
}

//---------------------------------------------------------------------------//
// end of Utils/utils/test/tstHyperslab_Vector.cc
//---------------------------------------------------------------------------//
