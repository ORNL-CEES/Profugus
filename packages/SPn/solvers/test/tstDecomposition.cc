//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/test/tstDecomposition.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 17 21:12:11 2014
 * \brief  Decomposition unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "gtest/utils_gtest.hh"

#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "../Decomposition.hh"


using namespace std;

//---------------------------------------------------------------------------//
// FIXTURE
//---------------------------------------------------------------------------//

class Decomposition_Test : public testing::Test
{
  protected:
    typedef profugus::Decomposition Decomposition;
    typedef Decomposition::Vec_Int  Vec_Int;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

  protected:
    int node, nodes;
};

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

TEST_F(Decomposition_Test, 1PE)
{
    if (nodes != 1)
        return;

    // use single processor constructor
    Decomposition d(16);

    EXPECT_EQ(16, d.num_local());
    EXPECT_EQ(16, d.num_global());

    EXPECT_EQ(0,  d.global(0));
    EXPECT_EQ(1,  d.global(1));
    EXPECT_EQ(2,  d.global(2));
    EXPECT_EQ(3,  d.global(3));
    EXPECT_EQ(4,  d.global(4));
    EXPECT_EQ(5,  d.global(5));
    EXPECT_EQ(6,  d.global(6));
    EXPECT_EQ(7,  d.global(7));
    EXPECT_EQ(8,  d.global(8));
    EXPECT_EQ(9,  d.global(9));
    EXPECT_EQ(10, d.global(10));
    EXPECT_EQ(11, d.global(11));
    EXPECT_EQ(12, d.global(12));
    EXPECT_EQ(13, d.global(13));
    EXPECT_EQ(14, d.global(14));
    EXPECT_EQ(15, d.global(15));
}

//---------------------------------------------------------------------------//
// decomposition:
//
//  |---|---|---|---|        |---|---|---|---|
//  | 12| 13| 14| 15|        | 6 | 7 | 6 | 7 |
//  |---|---|---|---|        |---|---|---|---|
//  | 8 | 9 | 10| 11|        | 4 | 5 | 4 | 5 |
//  |---|---|---|---|        |---|---|---|---|
//  | 4 | 5 | 6 | 7 |        | 2 | 3 | 2 | 3 |
//  |---|---|---|---|        |---|---|---|---|
//  | 0 | 1 | 2 | 3 |        | 0 | 1 | 0 | 1 |
//  |---|---|---|---|        |---|---|---|---|

TEST_F(Decomposition_Test, 2PE)
{
    if (nodes != 2)
        return;

    // make a map
    Decomposition::Vec_Int ltg(8);

    int i_off = 0;
    int j_off = 0;

    if (node == 1)
    {
        i_off = 2;
    }

    for (int j = 0; j < 4; j++)
    {
        for (int i = 0; i < 2; i++)
        {
            int local  = i + j * 2;
            int global = (i + i_off) + (j + j_off) * 4;

            ltg[local] = global;
        }
    }

    // make the decomposition
    Decomposition d(ltg);

    EXPECT_EQ(8, d.num_local());
    EXPECT_EQ(16, d.num_global());

    if (node == 0)
    {
        EXPECT_EQ(0, d.global(0));
        EXPECT_EQ(1, d.global(1));
        EXPECT_EQ(4, d.global(2));
        EXPECT_EQ(5, d.global(3));
        EXPECT_EQ(8, d.global(4));
        EXPECT_EQ(9, d.global(5));
        EXPECT_EQ(12, d.global(6));
        EXPECT_EQ(13, d.global(7));
    }

    else if (node == 1)
    {
        EXPECT_EQ(2, d.global(0));
        EXPECT_EQ(3, d.global(1));
        EXPECT_EQ(6, d.global(2));
        EXPECT_EQ(7, d.global(3));
        EXPECT_EQ(10, d.global(4));
        EXPECT_EQ(11, d.global(5));
        EXPECT_EQ(14, d.global(6));
        EXPECT_EQ(15, d.global(7));
    }

    // make an Epetra vector and try some communication
    Epetra_Vector x(d.map());
    {
        double *view;
        x.ExtractView(&view);
        fill(&view[0], &view[x.MyLength()], 1.0 + node);
    }

    // EXPORTING
    {
        // make target maps
        vector<int> exports;
        if (node == 0)
        {
            exports.push_back(2);
            exports.push_back(6);
            exports.push_back(10);
            exports.push_back(14);
        }
        else if (node == 1)
        {
            exports.push_back(1);
            exports.push_back(5);
            exports.push_back(9);
            exports.push_back(13);
        }

        // make target map for exporting
        Decomposition::Map target(-1, exports.size(), &exports[0], 0, d.comm());
        Epetra_Export exporter(d.map(), target);

        // results of export
        Epetra_Vector y(target);

        y.Export(x, exporter, Insert);

        for (int i = 0; i < y.MyLength(); i++)
        {
            if (node == 0)
            {
                EXPECT_EQ(2.0, y[i]);
            }
            else
            {
                EXPECT_EQ(1.0, y[i]);
            }
        }
    }

    // IMPORTING
    {
        // make target maps
        vector<int> imports;
        if (node == 0)
        {
            imports.push_back(2);
            imports.push_back(6);
            imports.push_back(10);
            imports.push_back(14);
        }
        else if (node == 1)
        {
            imports.push_back(1);
            imports.push_back(5);
            imports.push_back(9);
            imports.push_back(13);
        }

        // make target map for exporting
        Decomposition::Map target(-1, imports.size(), &imports[0], 0, d.comm());
        Epetra_Export importer(target, d.map());

        // results of export
        Epetra_Vector y(target);

        y.Import(x, importer, Insert);

        for (int i = 0; i < y.MyLength(); i++)
        {
            if (node == 0)
            {
                EXPECT_EQ(2.0, y[i]);
            }
            else
            {
                EXPECT_EQ(1.0, y[i]);
            }
        }
    }
}

//---------------------------------------------------------------------------//
// decomposition:
//
//  |---|---|---|---|        |---|---|---|---|
//  | 12| 13| 14| 15|        | 2 | 3 | 2 | 3 |
//  |---|---|---|---|     1  |---|---|---|---|
//  | 8 | 9 | 10| 11|        | 0 | 1 | 0 | 1 |
//  |---|---|---|---|        |---|---|---|---|
//  | 4 | 5 | 6 | 7 |        | 2 | 3 | 2 | 3 |
//  |---|---|---|---|     0  |---|---|---|---|
//  | 0 | 1 | 2 | 3 |        | 0 | 1 | 0 | 1 |
//  |---|---|---|---|        |---|---|---|---|

TEST_F(Decomposition_Test, 4PE)
{
    if (nodes != 4)
        return;

    // make a map
    Decomposition::Vec_Int ltg(4);

    int i_off = 0;
    int j_off = 0;

    if (node == 1)
    {
        i_off = 2;
    }
    else if (node == 2)
    {
        j_off = 2;
    }
    else if (node == 3)
    {
        i_off = 2;
        j_off = 2;
    }

    for (int j = 0; j < 2; j++)
    {
        for (int i = 0; i < 2; i++)
        {
            int local  = i + j * 2;
            int global = (i + i_off) + (j + j_off) * 4;

            ltg[local] = global;
        }
    }

    // make the decomposition
    Decomposition d(ltg);

    EXPECT_EQ(4, d.num_local());
    EXPECT_EQ(16, d.num_global());

    if (node == 0)
    {
        EXPECT_EQ(0, d.global(0));
        EXPECT_EQ(1, d.global(1));
        EXPECT_EQ(4, d.global(2));
        EXPECT_EQ(5, d.global(3));
    }
    else if (node == 1)
    {
        EXPECT_EQ(2, d.global(0));
        EXPECT_EQ(3, d.global(1));
        EXPECT_EQ(6, d.global(2));
        EXPECT_EQ(7, d.global(3));
    }
    else if (node == 2)
    {
        EXPECT_EQ(8, d.global(0));
        EXPECT_EQ(9, d.global(1));
        EXPECT_EQ(12, d.global(2));
        EXPECT_EQ(13, d.global(3));
    }
    else if (node == 3)
    {
        EXPECT_EQ(10, d.global(0));
        EXPECT_EQ(11, d.global(1));
        EXPECT_EQ(14, d.global(2));
        EXPECT_EQ(15, d.global(3));
    }
}

//---------------------------------------------------------------------------//
//                 end of tstDecomposition.cc
//---------------------------------------------------------------------------//
