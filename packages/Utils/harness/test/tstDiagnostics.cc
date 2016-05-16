//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/harness/test/tstDiagnostics.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 03 21:47:18 2013
 * \brief  Diagnostics unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "../Soft_Equivalence.hh"
#include "../Diagnostics.hh"

using namespace profugus;

//---------------------------------------------------------------------------//
// TEST FILTER
//---------------------------------------------------------------------------//

class DiagnosticsTest : public testing::Test
{
  protected:
    void SetUp()
    {
        // add an integer to global value
        Diagnostics::integers["B"] = 51;
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(DiagnosticsTest, ints)
{
    // add an integer quantity
    Diagnostics::integers["A"];
    EXPECT_EQ(0, Diagnostics::integers["A"]);

    // now update some values
    Diagnostics::integers["A"] = 21;
    Diagnostics::integers["A"]++;

    EXPECT_EQ(22, Diagnostics::integers["A"]);
}

//---------------------------------------------------------------------------//

TEST_F(DiagnosticsTest, floats)
{
    Diagnostics::integers["A"] = 22;

    // make a fraction entry
    Diagnostics::doubles["A_of_B"] = Diagnostics::integers["A"] /
                                     static_cast<double>(
                                         Diagnostics::integers["B"]);

    // check it
    EXPECT_DOUBLE_EQ(22.0/51.0, Diagnostics::doubles["A_of_B"]);

    // erase and check
    EXPECT_EQ(1u, Diagnostics::doubles.erase("A_of_B"));
    EXPECT_EQ(0u, Diagnostics::doubles.count("A_of_B"));

    Diagnostics::integers.erase("A");
}

//---------------------------------------------------------------------------//

TEST_F(DiagnosticsTest, vectors)
{
    Diagnostics::vec_integers["A"];
    Diagnostics::vec_integers["B"];

    EXPECT_TRUE(Diagnostics::vec_integers["A"].empty());
    EXPECT_TRUE(Diagnostics::vec_integers["B"].empty());

    Diagnostics::vec_integers["A"].resize(10);
    EXPECT_EQ(10u, Diagnostics::vec_integers["A"].size());

    Diagnostics::vec_doubles["B"].resize(2);
    Diagnostics::vec_doubles["B"][0] = 1.1;
    Diagnostics::vec_doubles["B"][1] = 2.4;

    EXPECT_EQ(2u, Diagnostics::vec_doubles["B"].size());

    std::vector<double> ref(2, 1.1);
    ref[1] = 2.4;
    EXPECT_VEC_SOFT_EQ(ref, Diagnostics::vec_doubles["B"]);

    Diagnostics::vec_doubles["B"].clear();
    EXPECT_TRUE(Diagnostics::vec_integers["B"].empty());

    EXPECT_EQ(51, Diagnostics::integers["B"]);
}

//---------------------------------------------------------------------------//

TEST_F(DiagnosticsTest, macro)
{
    int level[4];
    level[0] = 1;
    level[1] = 0;
    level[2] = 0;
    level[3] = 0;

#ifdef UTILS_DIAGNOSTICS_LEVEL_1
    level[1] = 1;
    level[0] = 0;
#endif

#ifdef UTILS_DIAGNOSTICS_LEVEL_2
    level[2] = 1;
    level[0] = 0;
#endif

#ifdef UTILS_DIAGNOSTICS_LEVEL_3
    level[3] = 1;
    level[0] = 0;
#endif

    DIAGNOSTICS_ONE(integers["L1"] = 1);
    DIAGNOSTICS_TWO(integers["L2"] = 1);
    DIAGNOSTICS_THREE(integers["L3"] = 1);

    if (level[0] == 1)
    {
        EXPECT_EQ(0u, Diagnostics::integers.count("L1"));
        EXPECT_EQ(0u, Diagnostics::integers.count("L2"));
        EXPECT_EQ(0u, Diagnostics::integers.count("L3"));
    }

    if (level[1] == 1)
    {
        EXPECT_EQ(1u, Diagnostics::integers.count("L1"));
    }

    if (level[2] == 1)
    {
        EXPECT_EQ(1u, Diagnostics::integers.count("L2"));
    }

    if (level[3] == 1)
    {
        EXPECT_EQ(1u, Diagnostics::integers.count("L3"));
    }
}

//---------------------------------------------------------------------------//
//                 end of tstDiagnostics.cc
//---------------------------------------------------------------------------//
