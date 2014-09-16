//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFission_Matrix_Processor.cc
 * \author Thomas M. Evans
 * \date   Fri Sep 12 11:09:04 2014
 * \brief  Fission_Matrix_Processor test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Fission_Matrix_Processor.hh"

#include "gtest/utils_gtest.hh"

#include <utility>
#include "comm/SpinLock.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Fission_Matrix_ProcessorTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Fission_Matrix_Processor Processor;
    typedef Processor::Ordered_Graph           Ordered_Graph;
    typedef Processor::Ordered_Matrix          Ordered_Matrix;
    typedef Processor::Sparse_Matrix           Sparse_Matrix;
    typedef Processor::Denominator             Denominator;
    typedef Processor::Idx                     Idx;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        local_denominator.resize(4);

        {
            Processor::Idx_Hash h(4);
            Sparse_Matrix m(local_matrix.bucket_count(), h);
            std::swap(m, local_matrix);
        }

        if (node == 0)
        {
            local_matrix[Idx(0, 0)] = 4.0;
            local_matrix[Idx(1, 0)] = 2.0;
            local_matrix[Idx(3, 0)] = 1.0;

            local_matrix[Idx(0, 1)] = 3.0;
            local_matrix[Idx(1, 1)] = 1.0;

            local_matrix[Idx(2, 2)] = 2.0;

            local_matrix[Idx(3, 3)] = 8.0;

            local_denominator[0] = 1.0;
            local_denominator[1] = 2.0;
            local_denominator[2] = 1.0;
            local_denominator[3] = 2.0;
        }
        if (node == 1)
        {
            local_matrix[Idx(0, 1)] = 4.0;
            local_matrix[Idx(1, 1)] = 2.0;
            local_matrix[Idx(2, 1)] = 1.0;

            local_matrix[Idx(0, 2)] = 3.0;
            local_matrix[Idx(1, 2)] = 1.0;
            local_matrix[Idx(2, 2)] = 7.0;

            local_matrix[Idx(3, 3)] = 6.0;

            local_denominator[1] = 1.0;
            local_denominator[2] = 2.0;
            local_denominator[3] = 3.0;
        }
        if (node == 2)
        {
            local_matrix[Idx(0, 3)] = 1.0;
            local_matrix[Idx(1, 3)] = 3.0;
            local_matrix[Idx(2, 3)] = 8.0;
            local_matrix[Idx(3, 3)] = 9.0;

            local_denominator[3] = 1.0;

        }
        if (node == 3)
        {
            local_matrix[Idx(1, 1)] = 6.0;
            local_matrix[Idx(2, 1)] = 4.0;
            local_matrix[Idx(3, 1)] = 1.0;

            local_denominator[1] = 1.0;
        }

        processor.build_matrix(local_matrix, local_denominator);
    }

  protected:
    // >>> DATA

    int node, nodes;

    // local fission matrix data
    Sparse_Matrix local_matrix;
    Denominator   local_denominator;

    // processor
    Processor processor;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_ProcessorTest, one_proc)
{
    if (nodes != 1)
        return;

    // test graph
    {
        const Ordered_Graph &g = processor.graph();
        EXPECT_EQ(3 + 2 + 1 + 1, g.size());

        Ordered_Graph ref = {Idx(0, 0),
                             Idx(0, 1),
                             Idx(1, 0),
                             Idx(1, 1),
                             Idx(2, 2),
                             Idx(3, 0),
                             Idx(3, 3)};

        EXPECT_EQ(ref, g);
    }

    // test matrix
    {
        const Ordered_Matrix &m = processor.matrix();
        EXPECT_EQ(3 + 2 + 1 + 1, m.size());

        Ordered_Matrix ref = {4.0,
                              1.5,
                              2.0,
                              0.5,
                              2.0,
                              1.0,
                              4.0};

        EXPECT_VEC_SOFTEQ(m, ref, 1.0e-14);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_ProcessorTest, two_proc)
{
    if (nodes != 2)
        return;

    // test graph
    {
        const Ordered_Graph &g = processor.graph();
        EXPECT_EQ(3 + 3 + 3 + 1, g.size());

        Ordered_Graph ref = {Idx(0, 0),
                             Idx(0, 1),
                             Idx(0, 2),
                             Idx(1, 0),
                             Idx(1, 1),
                             Idx(1, 2),
                             Idx(2, 1),
                             Idx(2, 2),
                             Idx(3, 0),
                             Idx(3, 3)};

        EXPECT_EQ(ref, g);
    }

    // test matrix
    {
        const Ordered_Matrix &m = processor.matrix();
        EXPECT_EQ(3 + 3 + 3 + 1, m.size());

        Ordered_Matrix ref = {4.0/1.0,
                              7.0/3.0,
                              3.0/3.0,
                              2.0/1.0,
                              3.0/3.0,
                              1.0/3.0,
                              1.0/3.0,
                              9.0/3.0,
                              1.0/1.0,
                              14.0/5.0};

        EXPECT_VEC_SOFTEQ(m, ref, 1.0e-14);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_ProcessorTest, four_proc)
{
    if (nodes != 4)
        return;

    // test matrix
    {
        const Ordered_Graph &g = processor.graph();
        EXPECT_EQ(3 + 4 + 3 + 4, g.size());

        Ordered_Graph ref = {Idx(0, 0),
                             Idx(0, 1),
                             Idx(0, 2),
                             Idx(0, 3),
                             Idx(1, 0),
                             Idx(1, 1),
                             Idx(1, 2),
                             Idx(1, 3),
                             Idx(2, 1),
                             Idx(2, 2),
                             Idx(2, 3),
                             Idx(3, 0),
                             Idx(3, 1),
                             Idx(3, 3)};

        EXPECT_EQ(ref, g);
    }

    // test matrix
    {
        const Ordered_Matrix &m = processor.matrix();
        EXPECT_EQ(3 + 4 + 3 + 4, m.size());

        // 0 = 1.0
        // 1 = 4.0
        // 2 = 3.0
        // 3 = 6.0

        Ordered_Matrix ref = {4.0/1.0,
                              7.0/4.0,
                              3.0/3.0,
                              1.0/6.0,
                              2.0/1.0,
                              9.0/4.0,
                              1.0/3.0,
                              3.0/6.0,
                              5.0/4.0,
                              9.0/3.0,
                              8.0/6.0,
                              1.0/1.0,
                              1.0/4.0,
                              23.0/6.0};

        EXPECT_VEC_SOFTEQ(m, ref, 1.0e-14);
    }

    // test reset
    processor.reset();
    const auto &g = processor.graph();
    const auto &m = processor.matrix();

    EXPECT_TRUE(g.empty());
    EXPECT_TRUE(m.empty());
}

//---------------------------------------------------------------------------//
/*
          0

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_1)
{
    if (nodes != 1)
        return;

    std::vector<int> pt = {Processor::NONE};
    std::vector<int> lc = {Processor::NONE};
    std::vector<int> rc = {Processor::NONE};
    std::vector<int> t  = {Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /
       1

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_2)
{
    if (nodes != 2)
        return;

    std::vector<int> pt = {Processor::NONE, 0};
    std::vector<int> lc = {1, Processor::NONE};
    std::vector<int> rc = {Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /   \
       1     2

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_3)
{
    if (nodes != 3)
        return;

    std::vector<int> pt = {Processor::NONE, 0, 0};
    std::vector<int> lc = {1, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /   \
       1     2
      /
     3

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_4)
{
    if (nodes != 4)
        return;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1};
    std::vector<int> lc = {1, 3, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, Processor::NONE, Processor::NONE,
                           Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /   \
       1     2
      / \
     3   4

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_5)
{
    if (nodes != 5)
        return;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1};
    std::vector<int> lc = {1, 3, Processor::NONE, Processor::NONE,
                           Processor::NONE};
    std::vector<int> rc = {2, 4, Processor::NONE, Processor::NONE,
                           Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /   \
       1     2
      / \   /
     3   4 5

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_6)
{
    if (nodes != 6)
        return;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1, 2};
    std::vector<int> lc = {1, 3, 5, Processor::NONE,
                           Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, 4, Processor::NONE, Processor::NONE,
                           Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::INTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /   \
       1     2
      / \   / \
     3   4 5   6

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_7)
{
    if (nodes != 7)
        return;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1, 2, 2};
    std::vector<int> lc = {1, 3, 5, Processor::NONE,
                           Processor::NONE, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, 4, 6, Processor::NONE, Processor::NONE,
                           Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::INTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//
/*
          0
        /   \
       1     2
      / \   / \
     3   4 5   6
    /
   7

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_8)
{
    if (nodes != 8)
        return;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1, 2, 2, 3};
    std::vector<int> lc = {1, 3, 5, 7, Processor::NONE,
                           Processor::NONE, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, 4, 6, Processor::NONE, Processor::NONE,
                           Processor::NONE, Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::INTERNAL, Processor::INTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  processor.node_type());
    EXPECT_EQ(pt[node], processor.parent());
    EXPECT_EQ(lc[node], processor.children()[0]);
    EXPECT_EQ(rc[node], processor.children()[1]);
}

//---------------------------------------------------------------------------//

TEST(Vector_Pair, continuous)
{
    typedef std::pair<int, int> Int_Pair;

    EXPECT_EQ(8, sizeof(Int_Pair));

    std::vector<Int_Pair> x(3);
    x[0] = Int_Pair(1, 2);
    x[1] = Int_Pair(3, 4);
    x[2] = Int_Pair(5, 6);

    int *start = &x[0].first;

    EXPECT_EQ(1, start[0]);
    EXPECT_EQ(2, start[1]);
    EXPECT_EQ(3, start[2]);
    EXPECT_EQ(4, start[3]);
    EXPECT_EQ(5, start[4]);
    EXPECT_EQ(6, start[5]);

    for (auto itr = &x[0].first; itr != (start + 6); ++itr)
    {
        *itr += 1;
    }

    EXPECT_EQ(Int_Pair(2, 3), x[0]);
    EXPECT_EQ(Int_Pair(4, 5), x[1]);
    EXPECT_EQ(Int_Pair(6, 7), x[2]);

    x.push_back(Int_Pair(6, 6));
    x.push_back(Int_Pair(8, 8));
    x.push_back(Int_Pair(0, 4));
    x.push_back(Int_Pair(8, 8));
    x.push_back(Int_Pair(8, 8));

    std::sort(x.begin(), x.end());

    auto insert = std::lower_bound(x.begin(), x.end(), Int_Pair(8, 8));

    EXPECT_EQ(8, x.size());
    EXPECT_EQ(5, insert - x.begin());
    EXPECT_EQ(Int_Pair(6, 7), x[4]);
    EXPECT_EQ(Int_Pair(8, 8), x[5]);
    EXPECT_EQ(Int_Pair(8, 8), x[6]);
    EXPECT_EQ(Int_Pair(8, 8), x[7]);

    auto done = std::lower_bound(x.begin(), x.end(), Int_Pair(9, 9));
    EXPECT_EQ(x.end(), done);
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Matrix_Processor.cc
//---------------------------------------------------------------------------//
