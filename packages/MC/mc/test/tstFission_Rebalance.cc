//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/test/tstFission_Rebalance.cc
 * \author Thomas M. Evans
 * \date   Mon May 05 12:19:09 2014
 * \brief  Fission_Rebalance unit-tests.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Fission_Rebalance.hh"

#include "gtest/utils_gtest.hh"
#include "geometry/RTK_Geometry.hh"

#include <memory>
#include <vector>
#include "utils/Definitions.hh"

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Fission_RebalanceTest : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::Core                          Geometry_t;
    typedef profugus::Fission_Rebalance<Geometry_t> Rebalance;
    typedef std::shared_ptr<Rebalance>              SP_Rebalance;
    typedef Rebalance::Fission_Site_t               Fission_Site_t;
    typedef Rebalance::Fission_Site_Container_t     Fission_Container_t;
    typedef Rebalance::Space_Vector                 Space_Vector;
    typedef Rebalance::Array_Bnds                   Array_Bnds;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        rebalance = std::make_shared<Rebalance>();
    }

    void setup(int N,
               int p0,
               int p1,
               int p2,
               int p3)
    {
        int offset = 0;

        if (node == 0)
        {
            offset = p0;

            // add 5 sites
            bank.resize(p1 - p0);
        }
        if (node == 1)
        {
            offset = p1;

            // add 5 sites
            bank.resize(p2 - p1);
        }
        if (node == 2)
        {
            offset = p2;

            // add 8 sites
            bank.resize(p3 - p2);
        }
        if (node == 3)
        {
            offset = p3;

            // add 7 sites
            bank.resize(N - p3);
        }

        for (int i = 0; i < bank.size(); ++i)
        {
            double pos = static_cast<double>(i + offset);
            bank[i].r  = Space_Vector(pos, pos, pos);
            bank[i].m  = i + offset;
        }
    }

    void check(int N)
    {
        std::vector<int> count(N, 0);
        int sum = 0;

        for (int n = 0; n < bank.size(); ++n)
        {
            ++sum;
            ++count[bank[n].m];

            EXPECT_EQ(bank[n].m, static_cast<int>(bank[n].r[0]));
            EXPECT_EQ(bank[n].m, static_cast<int>(bank[n].r[1]));
            EXPECT_EQ(bank[n].m, static_cast<int>(bank[n].r[2]));
        }

        profugus::global_sum(&count[0], N);
        profugus::global_sum(sum);

        EXPECT_EQ(N, sum);

        for (int n = 0; n < N; ++n)
        {
            EXPECT_EQ(1, count[n]);
        }
    }

  protected:
    // >>> Data that get re-initialized between tests

    SP_Rebalance rebalance;

    Fission_Container_t bank;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Fission_RebalanceTest, One_Iteration)
{
    setup(25, 0, 5, 10, 18);

    rebalance->rebalance(bank);
    const Array_Bnds &tb = rebalance->target_array_bnds();

    if (nodes == 1)
    {
        EXPECT_EQ(0, tb.first);
        EXPECT_EQ(0, tb.second);
        EXPECT_EQ(0, rebalance->num_iterations());
        EXPECT_EQ(5, bank.size());
        EXPECT_EQ(5, rebalance->num_fissions());
        EXPECT_EQ(5, rebalance->num_global_fissions());
    }

    if (nodes != 4) return;

    EXPECT_EQ(1, rebalance->num_iterations());

    EXPECT_EQ(25, rebalance->num_global_fissions());

    if (node == 0)
    {
        EXPECT_EQ(0, tb.first);
        EXPECT_EQ(6, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(7, bank.size());
        EXPECT_EQ(7, rebalance->num_fissions());
    }
    if (node == 1)
    {
        EXPECT_EQ(7, tb.first);
        EXPECT_EQ(12, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(6, bank.size());
        EXPECT_EQ(6, rebalance->num_fissions());
    }
    if (node == 2)
    {
        EXPECT_EQ(13, tb.first);
        EXPECT_EQ(18, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(6, bank.size());
        EXPECT_EQ(6, rebalance->num_fissions());
    }
    if (node == 3)
    {
        EXPECT_EQ(19, tb.first);
        EXPECT_EQ(24, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(0, rebalance->num_receives());

        EXPECT_EQ(6, bank.size());
        EXPECT_EQ(6, rebalance->num_fissions());
    }

    check(25);
}

//---------------------------------------------------------------------------//

TEST_F(Fission_RebalanceTest, Two_Iteration)
{
    if (nodes != 4) return;

    setup(43, 0, 6, 9, 37);

    rebalance->rebalance(bank);

    EXPECT_EQ(43, rebalance->num_global_fissions());

    const Array_Bnds &tb = rebalance->target_array_bnds();
    EXPECT_EQ(2, rebalance->num_iterations());

    if (node == 0)
    {
        EXPECT_EQ(0, tb.first);
        EXPECT_EQ(10, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(2, rebalance->num_receives());

        EXPECT_EQ(11, bank.size());
    }
    if (node == 1)
    {
        EXPECT_EQ(11, tb.first);
        EXPECT_EQ(21, tb.second);

        EXPECT_EQ(2, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(11, bank.size());
    }
    if (node == 2)
    {
        EXPECT_EQ(22, tb.first);
        EXPECT_EQ(32, tb.second);

        EXPECT_EQ(2, rebalance->num_sends());
        EXPECT_EQ(0, rebalance->num_receives());

        EXPECT_EQ(11, bank.size());
    }
    if (node == 3)
    {
        EXPECT_EQ(33, tb.first);
        EXPECT_EQ(42, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(10, bank.size());
    }

    check(43);
}

//---------------------------------------------------------------------------//

TEST_F(Fission_RebalanceTest, High_Left)
{
    if (nodes != 4) return;

    setup(100, 0, 100, 100, 100);

    rebalance->rebalance(bank);

    const Array_Bnds &tb = rebalance->target_array_bnds();
    EXPECT_EQ(3, rebalance->num_iterations());

    if (node == 0)
    {
        EXPECT_EQ(0, tb.first);
        EXPECT_EQ(24, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(0, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }
    if (node == 1)
    {
        EXPECT_EQ(25, tb.first);
        EXPECT_EQ(49, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }
    if (node == 2)
    {
        EXPECT_EQ(50, tb.first);
        EXPECT_EQ(74, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }
    if (node == 3)
    {
        EXPECT_EQ(75, tb.first);
        EXPECT_EQ(99, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }

    check(100);
}

//---------------------------------------------------------------------------//

TEST_F(Fission_RebalanceTest, High_Right)
{
    if (nodes != 4) return;

    setup(102, 0, 0, 0, 0);

    rebalance->rebalance(bank);

    const Array_Bnds &tb = rebalance->target_array_bnds();
    EXPECT_EQ(3, rebalance->num_iterations());

    if (node == 0)
    {
        EXPECT_EQ(0, tb.first);
        EXPECT_EQ(25, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(26, bank.size());
    }
    if (node == 1)
    {
        EXPECT_EQ(26, tb.first);
        EXPECT_EQ(51, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(26, bank.size());
    }
    if (node == 2)
    {
        EXPECT_EQ(52, tb.first);
        EXPECT_EQ(76, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }
    if (node == 3)
    {
        EXPECT_EQ(77, tb.first);
        EXPECT_EQ(101, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(0, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }

    check(102);
}

//---------------------------------------------------------------------------//

TEST_F(Fission_RebalanceTest, Multi_Rebalance)
{
    if (nodes != 4) return;

    // first rebalance
    setup(100, 0, 34, 51, 87);
    rebalance->rebalance(bank);
    check(100);

    EXPECT_EQ(1, rebalance->num_iterations());
    EXPECT_EQ(25, bank.size());

    // second rebalance
    setup(97, 0, 24, 49, 78);
    rebalance->rebalance(bank);
    const Array_Bnds &tb = rebalance->target_array_bnds();
    check(97);

    EXPECT_EQ(1, rebalance->num_iterations());

    if (node == 0)
    {
        EXPECT_EQ(0, tb.first);
        EXPECT_EQ(24, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(25, bank.size());
    }
    if (node == 1)
    {
        EXPECT_EQ(25, tb.first);
        EXPECT_EQ(48, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(0, rebalance->num_receives());

        EXPECT_EQ(24, bank.size());
    }
    if (node == 2)
    {
        EXPECT_EQ(49, tb.first);
        EXPECT_EQ(72, tb.second);

        EXPECT_EQ(1, rebalance->num_sends());
        EXPECT_EQ(0, rebalance->num_receives());

        EXPECT_EQ(24, bank.size());
    }
    if (node == 3)
    {
        EXPECT_EQ(73, tb.first);
        EXPECT_EQ(96, tb.second);

        EXPECT_EQ(0, rebalance->num_sends());
        EXPECT_EQ(1, rebalance->num_receives());

        EXPECT_EQ(24, bank.size());
    }

    // third rebalance
    setup(104, 0, 8, 19, 98);
    rebalance->rebalance(bank);
    check(104);

    EXPECT_EQ(2, rebalance->num_iterations());
    EXPECT_EQ(26, bank.size());
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Rebalance.cc
//---------------------------------------------------------------------------//
