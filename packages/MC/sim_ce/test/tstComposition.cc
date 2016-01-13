//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/sim_ce/test/tstComposition.cc
 * \author Thomas M Evans
 * \date   Mon Jan 11 15:24:34 2016
 * \brief  Composition class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Composition.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class CompositionTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    using Composition = profugus::Composition;
    using Vec_Zaids   = Composition::Vec_Zaids;
    using Vec_Dbl     = Composition::Vec_Dbl;

  protected:
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CompositionTest, construction)
{
    Composition a, b;

    // construct a
    {
        Vec_Zaids z = {1001, 8016};
        Vec_Dbl   N = {0.06684, 0.03342};

        a.add(z, N);
    }

    // construct b (with move)
    {
        Vec_Zaids z = {8016, 92234, 92235, 92236, 92238};
        Vec_Dbl   N = {0.046227, 7.49947e-06, 0.000865724, 3.96543e-06,
                       0.0222363};

        b.add(std::move(z), std::move(N));
    }

    EXPECT_EQ(2, a.num_nuclides());
    EXPECT_EQ(5, b.num_nuclides());

    Vec_Zaids rza = {1001, 8016};
    Vec_Dbl   rNa = {0.06684, 0.03342};
    Vec_Zaids rzb = {8016, 92234, 92235, 92236, 92238};
    Vec_Dbl   rNb = {0.046227, 7.49947e-06, 0.000865724, 3.96543e-06,
                     0.0222363};

    EXPECT_VEC_EQ(rza, a.zaids());
    EXPECT_VEC_EQ(rzb, b.zaids());

    EXPECT_VEC_EQ(rNa, a.number_dens());
    EXPECT_VEC_EQ(rNb, b.number_dens());

    for (int n = 0; n < 2; ++n)
    {
        EXPECT_EQ(rza[n], a.zaid(n));
        EXPECT_EQ(rNa[n], a.number_den(n));
    }
}

//---------------------------------------------------------------------------//
// end of MC/sim_ce/test/tstComposition.cc
//---------------------------------------------------------------------------//
