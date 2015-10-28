//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstBank.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 16:50:26 2014
 * \brief  Bank test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Bank.hh"

//---------------------------------------------------------------------------//
// TESTS (simple test harness)
//---------------------------------------------------------------------------//

class BankTest : public ::testing::Test
{
  protected:
    typedef profugus::Bank      Bank_t;
    typedef Bank_t::Particle_t  Particle_t;

    void SetUp()
    {
        // Initialization that are performed for each test
        // changes to these don't propagate between tests
        m_orig_p.set_wt(1.23);
        m_orig_p.set_matid(1);
    }

    // data that get re-initialized between tests
    Particle_t m_orig_p;
    Bank_t b;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(BankTest, empty)
{
    // empty bank
    EXPECT_TRUE(b.empty());
    EXPECT_EQ(0, b.size());
    EXPECT_EQ(0, b.num_particles());
    EXPECT_EQ(0, b.num_unique());
    EXPECT_EQ(0, b.next_count());
}

//---------------------------------------------------------------------------//

TEST_F(BankTest, one_particle)
{
    Particle_t p = m_orig_p;

    b.push(p);

    EXPECT_TRUE(!b.empty());
    EXPECT_EQ(1, b.size());
    EXPECT_EQ(1, b.num_particles());
    EXPECT_EQ(1, b.num_unique());

    // change our original particle
    p.set_wt(3.1415);

    // test the particle on top of the stack, see if it has the orig weight
    EXPECT_EQ(1.23, b.top().wt());

    // pop one particle
    auto popped = b.pop();

    EXPECT_EQ(1.23, popped.wt());

    // empty bank
    EXPECT_TRUE(b.empty());
    EXPECT_EQ(0, b.size());
    EXPECT_EQ(0, b.num_unique());
    EXPECT_EQ(0, b.num_particles());
}

//---------------------------------------------------------------------------//

TEST_F(BankTest, two_copies)
{
    // add two copies of one particle
    b.push(m_orig_p, 2u);

    EXPECT_TRUE(!b.empty());
    EXPECT_EQ(2, b.size());
    EXPECT_EQ(2, b.num_particles());
    EXPECT_EQ(1, b.num_unique());

    // pop two particles
    auto popped = b.pop();
    EXPECT_EQ(1.23, popped.wt());

    auto popped2 = b.pop();
    EXPECT_EQ(1.23, popped2.wt());

    // empty bank
    EXPECT_TRUE(b.empty());
    EXPECT_EQ(0, b.size());
    EXPECT_EQ(0, b.num_unique());
    EXPECT_EQ(0, b.num_particles());
}

//---------------------------------------------------------------------------//

TEST_F(BankTest, two_copies_of_two)
{
    Particle_t orig_p2;
    orig_p2.set_wt(1.23);
    orig_p2.set_matid(2);

    // add two copies of two particles
    b.push(m_orig_p, 2u);
    b.push(orig_p2, 2u);

    EXPECT_TRUE(!b.empty());
    EXPECT_EQ(4, b.size());
    EXPECT_EQ(4, b.num_particles());
    EXPECT_EQ(2, b.num_unique());

    // pop two particles
    EXPECT_EQ(2, b.pop().matid());
    EXPECT_EQ(2, b.pop().matid());
    EXPECT_EQ(1, b.pop().matid());
    EXPECT_EQ(1, b.pop().matid());

    // empty bank
    EXPECT_TRUE(b.empty());
    EXPECT_EQ(0, b.size());
    EXPECT_EQ(0, b.num_unique());
    EXPECT_EQ(0, b.num_particles());
}

//---------------------------------------------------------------------------//
//                 end of tstBank.cc
//---------------------------------------------------------------------------//
