//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstLG_Indexer.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 12 10:07:37 2014
 * \brief  LG_Indexer unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <vector>
#include <cmath>

#include "../LG_Indexer.hh"

using namespace std;

using profugus::LG_Indexer;

using def::I;
using def::J;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// decomposition:
//
//  |---|---|---|---|
//  | 12| 13| 14| 15|
//  |---|---|---|---|
//  | 8 | 9 | 10| 11|
//  |---|---|---|---|
//  | 4 | 5 | 6 | 7 |
//  |---|---|---|---|
//  | 0 | 1 | 2 | 3 |
//  |---|---|---|---|

TEST(LGIndexerOneProc, initializing)
{
    if (profugus::nodes() != 1)
    {
        SUCCEED() << "This test requires one processor";
        return;
    }

    vector<int> num_I(1, 4);
    vector<int> num_J(1, 4);

    LG_Indexer indexer(num_I, num_J);

    EXPECT_EQ(1, indexer.num_blocks(def::I));
    EXPECT_EQ(1, indexer.num_blocks(def::J));

    // the local and global indices should be the same
    EXPECT_EQ(0, indexer.l2g(0, 0, 0));
    EXPECT_EQ(1, indexer.l2g(1, 0, 0));
    EXPECT_EQ(2, indexer.l2g(2, 0, 0));
    EXPECT_EQ(3, indexer.l2g(3, 0, 0));
    EXPECT_EQ(4, indexer.l2g(0, 1, 0));
    EXPECT_EQ(5, indexer.l2g(1, 1, 0));
    EXPECT_EQ(6, indexer.l2g(2, 1, 0));
    EXPECT_EQ(7, indexer.l2g(3, 1, 0));
    EXPECT_EQ(8, indexer.l2g(0, 2, 0));
    EXPECT_EQ(9, indexer.l2g(1, 2, 0));
    EXPECT_EQ(10, indexer.l2g(2, 2, 0));
    EXPECT_EQ(11, indexer.l2g(3, 2, 0));
    EXPECT_EQ(12, indexer.l2g(0, 3, 0));
    EXPECT_EQ(13, indexer.l2g(1, 3, 0));
    EXPECT_EQ(14, indexer.l2g(2, 3, 0));
    EXPECT_EQ(15, indexer.l2g(3, 3, 0));

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            EXPECT_TRUE(indexer.convert_to_local(i, j) ==
                      LG_Indexer::IJ_Set(i, j));
            EXPECT_TRUE(indexer.convert_to_global(i, j) ==
                      LG_Indexer::IJ_Set(i, j));
        }
    }

    int i, j, k;
    indexer.l2l(0, i, j, k);
    indexer.l2l(0, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
    indexer.l2l(1, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
    indexer.l2l(2, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
    indexer.l2l(3, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);

    indexer.l2l(4, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
    indexer.l2l(5, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
    indexer.l2l(6, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
    indexer.l2l(7, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);

    indexer.l2l(8, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);
    indexer.l2l(9, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);
    indexer.l2l(10, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);
    indexer.l2l(11, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);

    indexer.l2l(12, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);
    indexer.l2l(13, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);
    indexer.l2l(14, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);
    indexer.l2l(15, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);

    EXPECT_EQ(4, indexer.num_cells(I));
    EXPECT_EQ(4, indexer.num_cells(J));

    EXPECT_EQ(4, indexer.num_cells(I));
    EXPECT_EQ(4, indexer.num_cells(J));
}

TEST(LGIndexerOneProc, global_to_manylocal)
{
    typedef LG_Indexer::IJ_Set IJ_Set;
    typedef LG_Indexer::Vec_IJ_Set Vec_IJ_Set;
    typedef std::vector<int> Vec_Int;

    if (profugus::nodes() != 1)
    {
        SUCCEED() << "This test requires one processor";
        return;
    }

    vector<int> num_I; num_I.push_back(4);
    vector<int> num_J; num_J.push_back(4);

    LG_Indexer indexer(num_I, num_J);
    Vec_Int    domains;
    Vec_IJ_Set lbegins;
    Vec_IJ_Set lends;

    // entire domain
    indexer.global_to_manylocal(IJ_Set(0,0), IJ_Set(4,4), 0,
            domains, lbegins, lends);

    EXPECT_EQ(1, domains.size());
    EXPECT_EQ(0, domains[0]);

    EXPECT_EQ(IJ_Set(0,0), lbegins[0]);
    EXPECT_EQ(IJ_Set(4,4), lends[0]);

    // part of the domain
    indexer.global_to_manylocal(IJ_Set(1,0), IJ_Set(3,2), 0,
            domains, lbegins, lends);

    EXPECT_EQ(1, domains.size());
    EXPECT_EQ(0, domains[0]);

    EXPECT_EQ(IJ_Set(1,0), lbegins[0]);
    EXPECT_EQ(IJ_Set(3,2), lends[0]);

}

//---------------------------------------------------------------------------//
// decomposition:
//
//  |---|---|---|---|     |---|---|---|---|
//  | 12| 13| 14| 15|     | 28| 29| 30| 31|
//  |---|---|---|---|     |---|---|---|---|
//  | 8 | 9 | 10| 11|     | 24| 25| 26| 27|
//  |---|---|---|---|     |---|---|---|---|
//  | 4 | 5 | 6 | 7 |     | 20| 21| 22| 23|
//  |---|---|---|---|     |---|---|---|---|
//  | 0 | 1 | 2 | 3 |     | 16| 17| 18| 19|
//  |---|---|---|---|     |---|---|---|---|
//
//
//  |---|---||---|---|     |---|---||---|---|
//  | 2 | 3 || 2 | 3 |     | 6 | 7 || 6 | 7 |
//  |---|---||---|---|     |---|---||---|---|    1
//  | 0 | 1 || 0 | 1 |     | 4 | 5 || 4 | 5 |
//  |=======||=======|     |=======||=======|
//  | 2 | 3 || 2 | 3 |     | 6 | 7 || 6 | 7 |
//  |---|---||---|---|     |---|---||---|---|    0
//  | 0 | 1 || 0 | 1 |     | 4 | 5 || 4 | 5 |
//  |---|---||---|---|     |---|---||---|---|
//
//      0        1             0        1

TEST(LGIndexer, FourProc)
{
    if (profugus::nodes() != 4)
    {
        SUCCEED() << "This test requires four processors";
        return;
    }

    int node = profugus::node();

    vector<int> num_I(2, 2);
    vector<int> num_J(2, 2);

    LG_Indexer indexer(num_I, num_J);

    EXPECT_EQ(2, indexer.num_blocks(def::I));
    EXPECT_EQ(2, indexer.num_blocks(def::J));

    EXPECT_EQ(4, indexer.num_global(def::I));
    EXPECT_EQ(4, indexer.num_global(def::J));

    EXPECT_EQ(4, indexer.num_blocks());
    EXPECT_EQ(0, indexer.domain(0, 0));
    EXPECT_EQ(1, indexer.domain(1, 0));
    EXPECT_EQ(2, indexer.domain(2, 0));
    EXPECT_EQ(3, indexer.domain(3, 0));

    if (node == 0)
    {
        EXPECT_EQ(0, indexer.l2g(0, 0, 0));
        EXPECT_EQ(1, indexer.l2g(1, 0, 0));
        EXPECT_EQ(4, indexer.l2g(0, 1, 0));
        EXPECT_EQ(5, indexer.l2g(1, 1, 0));
        EXPECT_EQ(16, indexer.l2g(0, 0, 1));
        EXPECT_EQ(17, indexer.l2g(1, 0, 1));
        EXPECT_EQ(20, indexer.l2g(0, 1, 1));
        EXPECT_EQ(21, indexer.l2g(1, 1, 1));

        EXPECT_EQ(2, indexer.num_cells(I));
        EXPECT_EQ(2, indexer.num_cells(J));

        EXPECT_TRUE(indexer.convert_to_global(0, 0) ==
                  LG_Indexer::IJ_Set(0, 0));
        EXPECT_TRUE(indexer.convert_to_global(1, 0) ==
                  LG_Indexer::IJ_Set(1, 0));
        EXPECT_TRUE(indexer.convert_to_global(0, 1) ==
                  LG_Indexer::IJ_Set(0, 1));
        EXPECT_TRUE(indexer.convert_to_global(1, 1) ==
                  LG_Indexer::IJ_Set(1, 1));
    }

    else if (node == 1)
    {
        EXPECT_EQ(2, indexer.l2g(0, 0, 0));
        EXPECT_EQ(3, indexer.l2g(1, 0, 0));
        EXPECT_EQ(6, indexer.l2g(0, 1, 0));
        EXPECT_EQ(7, indexer.l2g(1, 1, 0));
        EXPECT_EQ(18, indexer.l2g(0, 0, 1));
        EXPECT_EQ(19, indexer.l2g(1, 0, 1));
        EXPECT_EQ(22, indexer.l2g(0, 1, 1));
        EXPECT_EQ(23, indexer.l2g(1, 1, 1));

        EXPECT_EQ(2, indexer.num_cells(I));
        EXPECT_EQ(2, indexer.num_cells(J));

        EXPECT_TRUE(indexer.convert_to_global(0, 0) ==
                  LG_Indexer::IJ_Set(2, 0));
        EXPECT_TRUE(indexer.convert_to_global(1, 0) ==
                  LG_Indexer::IJ_Set(3, 0));
        EXPECT_TRUE(indexer.convert_to_global(0, 1) ==
                  LG_Indexer::IJ_Set(2, 1));
        EXPECT_TRUE(indexer.convert_to_global(1, 1) ==
                  LG_Indexer::IJ_Set(3, 1));
    }

    else if (node == 2)
    {
        EXPECT_EQ(8, indexer.l2g(0, 0, 0));
        EXPECT_EQ(9, indexer.l2g(1, 0, 0));
        EXPECT_EQ(12, indexer.l2g(0, 1, 0));
        EXPECT_EQ(13, indexer.l2g(1, 1, 0));
        EXPECT_EQ(24, indexer.l2g(0, 0, 1));
        EXPECT_EQ(25, indexer.l2g(1, 0, 1));
        EXPECT_EQ(28, indexer.l2g(0, 1, 1));
        EXPECT_EQ(29, indexer.l2g(1, 1, 1));

        EXPECT_EQ(2, indexer.num_cells(I));
        EXPECT_EQ(2, indexer.num_cells(J));

        EXPECT_TRUE(indexer.convert_to_global(0, 0) ==
                  LG_Indexer::IJ_Set(0, 2));
        EXPECT_TRUE(indexer.convert_to_global(1, 0) ==
                  LG_Indexer::IJ_Set(1, 2));
        EXPECT_TRUE(indexer.convert_to_global(0, 1) ==
                  LG_Indexer::IJ_Set(0, 3));
        EXPECT_TRUE(indexer.convert_to_global(1, 1) ==
                  LG_Indexer::IJ_Set(1, 3));
    }

    else if (node == 3)
    {
        EXPECT_EQ(10, indexer.l2g(0, 0, 0));
        EXPECT_EQ(11, indexer.l2g(1, 0, 0));
        EXPECT_EQ(14, indexer.l2g(0, 1, 0));
        EXPECT_EQ(15, indexer.l2g(1, 1, 0));
        EXPECT_EQ(26, indexer.l2g(0, 0, 1));
        EXPECT_EQ(27, indexer.l2g(1, 0, 1));
        EXPECT_EQ(30, indexer.l2g(0, 1, 1));
        EXPECT_EQ(31, indexer.l2g(1, 1, 1));

        EXPECT_EQ(2, indexer.num_cells(I));
        EXPECT_EQ(2, indexer.num_cells(J));

        EXPECT_TRUE(indexer.convert_to_global(0, 0) ==
                  LG_Indexer::IJ_Set(2, 2));
        EXPECT_TRUE(indexer.convert_to_global(1, 0) ==
                  LG_Indexer::IJ_Set(3, 2));
        EXPECT_TRUE(indexer.convert_to_global(0, 1) ==
                  LG_Indexer::IJ_Set(2, 3));
        EXPECT_TRUE(indexer.convert_to_global(1, 1) ==
                  LG_Indexer::IJ_Set(3, 3));
    }

    EXPECT_EQ(0, indexer.g2g(0, 0, 0));
    EXPECT_EQ(1, indexer.g2g(1, 0, 0));
    EXPECT_EQ(2, indexer.g2g(2, 0, 0));
    EXPECT_EQ(3, indexer.g2g(3, 0, 0));
    EXPECT_EQ(4, indexer.g2g(0, 1, 0));
    EXPECT_EQ(5, indexer.g2g(1, 1, 0));
    EXPECT_EQ(6, indexer.g2g(2, 1, 0));
    EXPECT_EQ(7, indexer.g2g(3, 1, 0));
    EXPECT_EQ(8, indexer.g2g(0, 2, 0));
    EXPECT_EQ(9, indexer.g2g(1, 2, 0));
    EXPECT_EQ(10, indexer.g2g(2, 2, 0));
    EXPECT_EQ(11, indexer.g2g(3, 2, 0));
    EXPECT_EQ(12, indexer.g2g(0, 3, 0));
    EXPECT_EQ(13, indexer.g2g(1, 3, 0));
    EXPECT_EQ(14, indexer.g2g(2, 3, 0));
    EXPECT_EQ(15, indexer.g2g(3, 3, 0));
    EXPECT_EQ(16, indexer.g2g(0, 0, 1));
    EXPECT_EQ(17, indexer.g2g(1, 0, 1));
    EXPECT_EQ(18, indexer.g2g(2, 0, 1));
    EXPECT_EQ(19, indexer.g2g(3, 0, 1));
    EXPECT_EQ(20, indexer.g2g(0, 1, 1));
    EXPECT_EQ(21, indexer.g2g(1, 1, 1));
    EXPECT_EQ(22, indexer.g2g(2, 1, 1));
    EXPECT_EQ(23, indexer.g2g(3, 1, 1));
    EXPECT_EQ(24, indexer.g2g(0, 2, 1));
    EXPECT_EQ(25, indexer.g2g(1, 2, 1));
    EXPECT_EQ(26, indexer.g2g(2, 2, 1));
    EXPECT_EQ(27, indexer.g2g(3, 2, 1));
    EXPECT_EQ(28, indexer.g2g(0, 3, 1));
    EXPECT_EQ(29, indexer.g2g(1, 3, 1));
    EXPECT_EQ(30, indexer.g2g(2, 3, 1));
    EXPECT_EQ(31, indexer.g2g(3, 3, 1));
}

TEST(LGIndexerFourProc, global_to_manylocal)
{
    typedef LG_Indexer::IJ_Set IJ_Set;
    typedef LG_Indexer::Vec_IJ_Set Vec_IJ_Set;
    typedef std::vector<int> Vec_Int;

    if (profugus::nodes() != 4)
    {
        SUCCEED() << "This test requires four processors";
        return;
    }

    int node = profugus::node();

    Vec_Int num_I(2, 2);
    Vec_Int num_J(2, 2);

    LG_Indexer indexer(num_I, num_J);
    Vec_Int    domains;
    Vec_IJ_Set lbegins;
    Vec_IJ_Set lends;

    indexer.global_to_manylocal(IJ_Set(0,0), IJ_Set(4,4), 0,
            domains, lbegins, lends);
    ASSERT_EQ(4, domains.size());
    EXPECT_EQ(0, domains[0]);
    EXPECT_EQ(1, domains[1]);
    EXPECT_EQ(2, domains[2]);
    EXPECT_EQ(3, domains[3]);

    EXPECT_EQ(IJ_Set(0,0), lbegins[0]);
    EXPECT_EQ(IJ_Set(2,2), lends[0]);
    EXPECT_EQ(IJ_Set(0,0), lbegins[3]);
    EXPECT_EQ(IJ_Set(2,2), lends[3]);

    indexer.global_to_manylocal(IJ_Set(1,1), IJ_Set(3,3), 0,
            domains, lbegins, lends);
    ASSERT_EQ(4, domains.size());

    EXPECT_EQ(IJ_Set(1,1), lbegins[0]);
    EXPECT_EQ(IJ_Set(2,2), lends[0]);
    EXPECT_EQ(IJ_Set(0,0), lbegins[3]);
    EXPECT_EQ(IJ_Set(1,1), lends[3]);
}

//---------------------------------------------------------------------------//

// decomposition:
//
//  |---|---|---|---|     |---|---|---|---|
//  | 12| 13| 14| 15|     | 28| 29| 30| 31|
//  |---|---|---|---|     |---|---|---|---|
//  | 8 | 9 | 10| 11|     | 24| 25| 26| 27|
//  |---|---|---|---|     |---|---|---|---|
//  | 4 | 5 | 6 | 7 |     | 20| 21| 22| 23|
//  |---|---|---|---|     |---|---|---|---|
//  | 0 | 1 | 2 | 3 |     | 16| 17| 18| 19|
//  |---|---|---|---|     |---|---|---|---|
//
//
//  |---|---|---||---|    |---|---|---||---|
//  | 3 | 4 | 5 || 1 |    | 9 | 10| 11|| 3 |
//  |---|---|---||---|    |---|---|---||---|    1
//  | 0 | 1 | 2 || 0 |    | 6 | 7 | 8 || 2 |
//  |===========||===|    |===========||===|
//  | 3 | 4 | 5 || 1 |    | 9 | 10| 11|| 3 |
//  |---|---|---||---|    |---|---|---||---|    0
//  | 0 | 1 | 2 || 0 |    | 6 | 7 | 8 || 2 |
//  |---|---|---||---|    |---|---|---||---|
//
//        0        1            0        1

TEST(LGIndexerFourProc, set_to_domain)
{
    if (profugus::nodes() != 4)
    {
        SUCCEED() << "This test requires one processor";
        return;
    }

    vector<int> num_I(2);
    vector<int> num_J(2, 2);
    num_I[0] = 3;
    num_I[1] = 1;

    LG_Indexer indexer(num_I, num_J);

    EXPECT_EQ(2, indexer.num_blocks(def::I));
    EXPECT_EQ(2, indexer.num_blocks(def::J));

    indexer.set_to_domain(0);

    EXPECT_EQ(4, indexer.num_global(def::I));
    EXPECT_EQ(4, indexer.num_global(def::J));

    EXPECT_EQ(0, indexer.offset(def::I));
    EXPECT_EQ(0, indexer.offset(def::J));

    EXPECT_EQ(0, indexer.current_domain());

    EXPECT_EQ(3, indexer.num_cells(I));
    EXPECT_EQ(2, indexer.num_cells(J));

    EXPECT_EQ(0, indexer.l2g(0, 0, 0));
    EXPECT_EQ(1, indexer.l2g(1, 0, 0));
    EXPECT_EQ(2, indexer.l2g(2, 0, 0));
    EXPECT_EQ(4, indexer.l2g(0, 1, 0));
    EXPECT_EQ(5, indexer.l2g(1, 1, 0));
    EXPECT_EQ(6, indexer.l2g(2, 1, 0));
    EXPECT_EQ(16, indexer.l2g(0, 0, 1));
    EXPECT_EQ(17, indexer.l2g(1, 0, 1));
    EXPECT_EQ(18, indexer.l2g(2, 0, 1));
    EXPECT_EQ(20, indexer.l2g(0, 1, 1));
    EXPECT_EQ(21, indexer.l2g(1, 1, 1));
    EXPECT_EQ(22, indexer.l2g(2, 1, 1));

    EXPECT_EQ(0, indexer.l2l(0, 0, 0));
    EXPECT_EQ(1, indexer.l2l(1, 0, 0));
    EXPECT_EQ(2, indexer.l2l(2, 0, 0));
    EXPECT_EQ(3, indexer.l2l(0, 1, 0));
    EXPECT_EQ(4, indexer.l2l(1, 1, 0));
    EXPECT_EQ(5, indexer.l2l(2, 1, 0));
    EXPECT_EQ(6, indexer.l2l(0, 0, 1));
    EXPECT_EQ(7, indexer.l2l(1, 0, 1));
    EXPECT_EQ(8, indexer.l2l(2, 0, 1));
    EXPECT_EQ(9, indexer.l2l(0, 1, 1));
    EXPECT_EQ(10, indexer.l2l(1, 1, 1));
    EXPECT_EQ(11, indexer.l2l(2, 1, 1));

    indexer.set_to_domain(3);

    EXPECT_EQ(3, indexer.current_domain());

    EXPECT_EQ(4, indexer.num_global(def::I));
    EXPECT_EQ(4, indexer.num_global(def::J));

    EXPECT_EQ(3, indexer.offset(def::I));
    EXPECT_EQ(2, indexer.offset(def::J));

    EXPECT_EQ(1, indexer.num_cells(I));
    EXPECT_EQ(2, indexer.num_cells(J));

    EXPECT_EQ(11, indexer.l2g(0, 0, 0));
    EXPECT_EQ(15, indexer.l2g(0, 1, 0));
    EXPECT_EQ(27, indexer.l2g(0, 0, 1));
    EXPECT_EQ(31, indexer.l2g(0, 1, 1));

    EXPECT_EQ(0, indexer.l2l(0, 0, 0));
    EXPECT_EQ(1, indexer.l2l(0, 1, 0));
    EXPECT_EQ(2, indexer.l2l(0, 0, 1));
    EXPECT_EQ(3, indexer.l2l(0, 1, 1));
}

//---------------------------------------------------------------------------//

TEST(LGIndexerFourProc, l2g)
{
    if (profugus::nodes() != 4)
    {
        SUCCEED() << "This test requires four processors";
        return;
    }
    int node = profugus::node();

    vector<int> num_I(2);
    vector<int> num_J(2, 2);
    num_I[0] = 3;
    num_I[1] = 1;

    LG_Indexer indexer(num_I, num_J);

    int i, j, k;

    if (node == 0)
    {
        indexer.l2l(0, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(1, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(2, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(3, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(4, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(5, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(6, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(7, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(8, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(9, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
        indexer.l2l(10, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
        indexer.l2l(11, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);

    }
    else if (node == 1)
    {
        indexer.l2l(0, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(1, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(2, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(3, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);

    }
    else if (node == 2)
    {
        indexer.l2l(0, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(1, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(2, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(3, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(4, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(5, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(6, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(7, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(8, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(9, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
        indexer.l2l(10, i, j, k);
        EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
        indexer.l2l(11, i, j, k);
        EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);


    }
    else if (node == 3)
    {
        indexer.l2l(0, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
        indexer.l2l(1, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
        indexer.l2l(2, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
        indexer.l2l(3, i, j, k);
        EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
    }

    indexer.g2g(0, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
    indexer.g2g(1, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
    indexer.g2g(2, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);
    indexer.g2g(3, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(0, j); EXPECT_EQ(0, k);

    indexer.g2g(4, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
    indexer.g2g(5, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
    indexer.g2g(6, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);
    indexer.g2g(7, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(1, j); EXPECT_EQ(0, k);

    indexer.g2g(8, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);
    indexer.g2g(9, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);
    indexer.g2g(10, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);
    indexer.g2g(11, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(2, j); EXPECT_EQ(0, k);

    indexer.g2g(12, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);
    indexer.g2g(13, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);
    indexer.g2g(14, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);
    indexer.g2g(15, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(3, j); EXPECT_EQ(0, k);

    indexer.g2g(16, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
    indexer.g2g(17, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
    indexer.g2g(18, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);
    indexer.g2g(19, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(0, j); EXPECT_EQ(1, k);

    indexer.g2g(20, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
    indexer.g2g(21, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
    indexer.g2g(22, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);
    indexer.g2g(23, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(1, j); EXPECT_EQ(1, k);

    indexer.g2g(24, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(2, j); EXPECT_EQ(1, k);
    indexer.g2g(25, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(2, j); EXPECT_EQ(1, k);
    indexer.g2g(26, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(2, j); EXPECT_EQ(1, k);
    indexer.g2g(27, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(2, j); EXPECT_EQ(1, k);

    indexer.g2g(28, i, j, k);
    EXPECT_EQ(0, i); EXPECT_EQ(3, j); EXPECT_EQ(1, k);
    indexer.g2g(29, i, j, k);
    EXPECT_EQ(1, i); EXPECT_EQ(3, j); EXPECT_EQ(1, k);
    indexer.g2g(30, i, j, k);
    EXPECT_EQ(2, i); EXPECT_EQ(3, j); EXPECT_EQ(1, k);
    indexer.g2g(31, i, j, k);
    EXPECT_EQ(3, i); EXPECT_EQ(3, j); EXPECT_EQ(1, k);
}

//---------------------------------------------------------------------------//

TEST(LGIndexerFourProc, convert_to_local)
{
    if (profugus::nodes() != 4)
    {
        SUCCEED() << "This test requires four processors";
        return;
    }
    int node = profugus::node();

    vector<int> num_I(2);
    vector<int> num_J(2, 2);
    num_I[0] = 3;
    num_I[1] = 1;

    LG_Indexer indexer(num_I, num_J);

    if (node == 0)
    {
        EXPECT_TRUE(indexer.convert_to_local(0, 0) ==
                  LG_Indexer::IJ_Set(0, 0));
        EXPECT_TRUE(indexer.convert_to_local(1, 0) ==
                  LG_Indexer::IJ_Set(1, 0));
        EXPECT_TRUE(indexer.convert_to_local(2, 0) ==
                  LG_Indexer::IJ_Set(2, 0));
        EXPECT_TRUE(indexer.convert_to_local(3, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 1) ==
                  LG_Indexer::IJ_Set(0, 1));
        EXPECT_TRUE(indexer.convert_to_local(1, 1) ==
                  LG_Indexer::IJ_Set(1, 1));
        EXPECT_TRUE(indexer.convert_to_local(2, 1) ==
                  LG_Indexer::IJ_Set(2, 1));
        EXPECT_TRUE(indexer.convert_to_local(3, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
    }

    else if (node == 1)
    {
        EXPECT_TRUE(indexer.convert_to_local(0, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 0) ==
                  LG_Indexer::IJ_Set(0, 0));

        EXPECT_TRUE(indexer.convert_to_local(0, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 1) ==
                  LG_Indexer::IJ_Set(0, 1));

        EXPECT_TRUE(indexer.convert_to_local(0, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
    }

    else if (node == 2)
    {
        EXPECT_TRUE(indexer.convert_to_local(0, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 2) ==
                  LG_Indexer::IJ_Set(0, 0));
        EXPECT_TRUE(indexer.convert_to_local(1, 2) ==
                  LG_Indexer::IJ_Set(1, 0));
        EXPECT_TRUE(indexer.convert_to_local(2, 2) ==
                  LG_Indexer::IJ_Set(2, 0));
        EXPECT_TRUE(indexer.convert_to_local(3, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 3) ==
                  LG_Indexer::IJ_Set(0, 1));
        EXPECT_TRUE(indexer.convert_to_local(1, 3) ==
                  LG_Indexer::IJ_Set(1, 1));
        EXPECT_TRUE(indexer.convert_to_local(2, 3) ==
                  LG_Indexer::IJ_Set(2, 1));
        EXPECT_TRUE(indexer.convert_to_local(3, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
    }

    else if (node == 3)
    {
        EXPECT_TRUE(indexer.convert_to_local(0, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 0) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 1) ==
                  LG_Indexer::IJ_Set(-1, -1));

        EXPECT_TRUE(indexer.convert_to_local(0, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 2) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 2) ==
                  LG_Indexer::IJ_Set(0, 0));

        EXPECT_TRUE(indexer.convert_to_local(0, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(1, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(2, 3) ==
                  LG_Indexer::IJ_Set(-1, -1));
        EXPECT_TRUE(indexer.convert_to_local(3, 3) ==
                  LG_Indexer::IJ_Set(0, 1));
    }
}

//---------------------------------------------------------------------------//
// decomposition:
//
//  |---|---|---|---|     |---|---|---|---|
//  | 12| 13| 14| 15|     | 28| 29| 30| 31|
//  |---|---|---|---|     |---|---|---|---|
//  | 8 | 9 | 10| 11|     | 24| 25| 26| 27|
//  |---|---|---|---|     |---|---|---|---|
//  | 4 | 5 | 6 | 7 |     | 20| 21| 22| 23|
//  |---|---|---|---|     |---|---|---|---|
//  | 0 | 1 | 2 | 3 |     | 16| 17| 18| 19|
//  |---|---|---|---|     |---|---|---|---|
//
//  BLOCK DECOMPOSITION
//
//  |---|---||---|---|     |---|---||---|---|
//  | 6 | 7 || 6 | 7 |     | 14| 15|| 14| 15|
//  |---|---||---|---|     |---|---||---|---|
//  | 4 | 5 || 4 | 5 |     | 12| 13|| 12| 13|
//  |---|---||---|---|     |---|---||---|---|  0
//  | 2 | 3 || 2 | 3 |     | 10| 11|| 10| 11|
//  |---|---||---|---|     |---|---||---|---|
//  | 0 | 1 || 0 | 1 |     | 8 | 9 || 8 | 9 |
//  |---|---||---|---|     |---|---||---|---|
//
//      0        1             0        1
//
//  SET DECOMPOSITION
//
//  |-------||-------|     |-------||-------|
//  |       ||       |     |       ||       |
//  |       ||       |     |       ||       |
//  |       ||       |     |       ||       |
//  |   0   ||   1   |     |   0   ||   1   |
//  |  [0]  ||  [1]  |     |  [2]  ||  [3]  |
//  |       ||       |     |       ||       |
//  |       ||       |     |       ||       |
//  |-------||-------|     |-------||-------|
//
//         set 0                  set 1

TEST(LGIndexerFourProc, sets)
{
    if (profugus::nodes() != 4)
    {
        SUCCEED() << "This test requires four processors";
        return;
    }
    int node = profugus::node();

    vector<int> num_I(2, 2);
    vector<int> num_J(1, 4);
    LG_Indexer indexer(num_I, num_J, 2);

    EXPECT_EQ(2, indexer.num_sets());
    EXPECT_EQ(2, indexer.num_blocks());
    EXPECT_EQ(2, indexer.num_blocks(I));
    EXPECT_EQ(1, indexer.num_blocks(J));
    EXPECT_EQ(node, indexer.current_domain());

    EXPECT_EQ(0, indexer.domain(0, 0));
    EXPECT_EQ(1, indexer.domain(1, 0));
    EXPECT_EQ(2, indexer.domain(0, 1));
    EXPECT_EQ(3, indexer.domain(1, 1));

    if (node == 0)
    {
        EXPECT_EQ(0, indexer.set());
        EXPECT_EQ(0, indexer.block());

        EXPECT_EQ(0, indexer.l2g(0,0,0));
        EXPECT_EQ(1, indexer.l2g(1,0,0));
        EXPECT_EQ(4, indexer.l2g(0,1,0));
        EXPECT_EQ(5, indexer.l2g(1,1,0));
        EXPECT_EQ(8, indexer.l2g(0,2,0));
        EXPECT_EQ(9, indexer.l2g(1,2,0));
        EXPECT_EQ(12, indexer.l2g(0,3,0));
        EXPECT_EQ(13, indexer.l2g(1,3,0));
        EXPECT_EQ(16, indexer.l2g(0,0,1));
        EXPECT_EQ(17, indexer.l2g(1,0,1));
        EXPECT_EQ(20, indexer.l2g(0,1,1));
        EXPECT_EQ(21, indexer.l2g(1,1,1));
        EXPECT_EQ(24, indexer.l2g(0,2,1));
        EXPECT_EQ(25, indexer.l2g(1,2,1));
        EXPECT_EQ(28, indexer.l2g(0,3,1));
        EXPECT_EQ(29, indexer.l2g(1,3,1));

        EXPECT_EQ(0, indexer.offset(I));
        EXPECT_EQ(0, indexer.offset(J));
    }
    else if (node == 1)
    {
        EXPECT_EQ(0, indexer.set());
        EXPECT_EQ(1, indexer.block());

        EXPECT_EQ(2, indexer.l2g(0,0,0));
        EXPECT_EQ(3, indexer.l2g(1,0,0));
        EXPECT_EQ(6, indexer.l2g(0,1,0));
        EXPECT_EQ(7, indexer.l2g(1,1,0));
        EXPECT_EQ(10, indexer.l2g(0,2,0));
        EXPECT_EQ(11, indexer.l2g(1,2,0));
        EXPECT_EQ(14, indexer.l2g(0,3,0));
        EXPECT_EQ(15, indexer.l2g(1,3,0));
        EXPECT_EQ(18, indexer.l2g(0,0,1));
        EXPECT_EQ(19, indexer.l2g(1,0,1));
        EXPECT_EQ(22, indexer.l2g(0,1,1));
        EXPECT_EQ(23, indexer.l2g(1,1,1));
        EXPECT_EQ(26, indexer.l2g(0,2,1));
        EXPECT_EQ(27, indexer.l2g(1,2,1));
        EXPECT_EQ(30, indexer.l2g(0,3,1));
        EXPECT_EQ(31, indexer.l2g(1,3,1));

        EXPECT_EQ(2, indexer.offset(I));
        EXPECT_EQ(0, indexer.offset(J));
    }
    else if (node == 2)
    {
        EXPECT_EQ(1, indexer.set());
        EXPECT_EQ(0, indexer.block());

        EXPECT_EQ(0, indexer.l2g(0,0,0));
        EXPECT_EQ(1, indexer.l2g(1,0,0));
        EXPECT_EQ(4, indexer.l2g(0,1,0));
        EXPECT_EQ(5, indexer.l2g(1,1,0));
        EXPECT_EQ(8, indexer.l2g(0,2,0));
        EXPECT_EQ(9, indexer.l2g(1,2,0));
        EXPECT_EQ(12, indexer.l2g(0,3,0));
        EXPECT_EQ(13, indexer.l2g(1,3,0));
        EXPECT_EQ(16, indexer.l2g(0,0,1));
        EXPECT_EQ(17, indexer.l2g(1,0,1));
        EXPECT_EQ(20, indexer.l2g(0,1,1));
        EXPECT_EQ(21, indexer.l2g(1,1,1));
        EXPECT_EQ(24, indexer.l2g(0,2,1));
        EXPECT_EQ(25, indexer.l2g(1,2,1));
        EXPECT_EQ(28, indexer.l2g(0,3,1));
        EXPECT_EQ(29, indexer.l2g(1,3,1));

        EXPECT_EQ(0, indexer.offset(I));
        EXPECT_EQ(0, indexer.offset(J));
    }
    else if (node == 3)
    {
        EXPECT_EQ(1, indexer.set());
        EXPECT_EQ(1, indexer.block());

        EXPECT_EQ(2, indexer.l2g(0,0,0));
        EXPECT_EQ(3, indexer.l2g(1,0,0));
        EXPECT_EQ(6, indexer.l2g(0,1,0));
        EXPECT_EQ(7, indexer.l2g(1,1,0));
        EXPECT_EQ(10, indexer.l2g(0,2,0));
        EXPECT_EQ(11, indexer.l2g(1,2,0));
        EXPECT_EQ(14, indexer.l2g(0,3,0));
        EXPECT_EQ(15, indexer.l2g(1,3,0));
        EXPECT_EQ(18, indexer.l2g(0,0,1));
        EXPECT_EQ(19, indexer.l2g(1,0,1));
        EXPECT_EQ(22, indexer.l2g(0,1,1));
        EXPECT_EQ(23, indexer.l2g(1,1,1));
        EXPECT_EQ(26, indexer.l2g(0,2,1));
        EXPECT_EQ(27, indexer.l2g(1,2,1));
        EXPECT_EQ(30, indexer.l2g(0,3,1));
        EXPECT_EQ(31, indexer.l2g(1,3,1));

        EXPECT_EQ(2, indexer.offset(I));
        EXPECT_EQ(0, indexer.offset(J));
    }

    // on all processors change to domain 3
    indexer.set_to_domain(3);
    EXPECT_EQ(3, indexer.current_domain());
    EXPECT_EQ(1, indexer.set());
    EXPECT_EQ(1, indexer.block());

    EXPECT_EQ(2, indexer.l2g(0,0,0));
    EXPECT_EQ(3, indexer.l2g(1,0,0));
    EXPECT_EQ(6, indexer.l2g(0,1,0));
    EXPECT_EQ(7, indexer.l2g(1,1,0));
    EXPECT_EQ(10, indexer.l2g(0,2,0));
    EXPECT_EQ(11, indexer.l2g(1,2,0));
    EXPECT_EQ(14, indexer.l2g(0,3,0));
    EXPECT_EQ(15, indexer.l2g(1,3,0));
    EXPECT_EQ(18, indexer.l2g(0,0,1));
    EXPECT_EQ(19, indexer.l2g(1,0,1));
    EXPECT_EQ(22, indexer.l2g(0,1,1));
    EXPECT_EQ(23, indexer.l2g(1,1,1));
    EXPECT_EQ(26, indexer.l2g(0,2,1));
    EXPECT_EQ(27, indexer.l2g(1,2,1));
    EXPECT_EQ(30, indexer.l2g(0,3,1));
    EXPECT_EQ(31, indexer.l2g(1,3,1));

    EXPECT_EQ(2, indexer.offset(I));
    EXPECT_EQ(0, indexer.offset(J));
}

//---------------------------------------------------------------------------//
//                 end of tstLG_Indexer.cc
//---------------------------------------------------------------------------//
