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

#include "comm/SpinLock.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Fission_Matrix_ProcessorTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Fission_Matrix_Processor Processor;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

  protected:
    // >>> DATA

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*
          0

 */
TEST_F(Fission_Matrix_ProcessorTest, setup_1)
{
    if (nodes != 1)
        return;

    Processor p;

    std::vector<int> pt = {Processor::NONE};
    std::vector<int> lc = {Processor::NONE};
    std::vector<int> rc = {Processor::NONE};
    std::vector<int> t  = {Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0};
    std::vector<int> lc = {1, Processor::NONE};
    std::vector<int> rc = {Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0, 0};
    std::vector<int> lc = {1, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1};
    std::vector<int> lc = {1, 3, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, Processor::NONE, Processor::NONE,
                           Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1};
    std::vector<int> lc = {1, 3, Processor::NONE, Processor::NONE,
                           Processor::NONE};
    std::vector<int> rc = {2, 4, Processor::NONE, Processor::NONE,
                           Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1, 2};
    std::vector<int> lc = {1, 3, 5, Processor::NONE,
                           Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, 4, Processor::NONE, Processor::NONE,
                           Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::INTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1, 2, 2};
    std::vector<int> lc = {1, 3, 5, Processor::NONE,
                           Processor::NONE, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, 4, 6, Processor::NONE, Processor::NONE,
                           Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::INTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
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

    Processor p;

    std::vector<int> pt = {Processor::NONE, 0, 0, 1, 1, 2, 2, 3};
    std::vector<int> lc = {1, 3, 5, 7, Processor::NONE,
                           Processor::NONE, Processor::NONE, Processor::NONE};
    std::vector<int> rc = {2, 4, 6, Processor::NONE, Processor::NONE,
                           Processor::NONE, Processor::NONE, Processor::NONE};
    std::vector<int> t  = {Processor::INTERNAL, Processor::INTERNAL,
                           Processor::INTERNAL, Processor::INTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL,
                           Processor::EXTERNAL, Processor::EXTERNAL};

    EXPECT_EQ(t[node],  p.node_type());
    EXPECT_EQ(pt[node], p.parent());
    EXPECT_EQ(lc[node], p.children()[0]);
    EXPECT_EQ(rc[node], p.children()[1]);
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Matrix_Processor.cc
//---------------------------------------------------------------------------//
